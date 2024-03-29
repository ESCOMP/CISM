!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_sparse.F90 - part of the Community Ice Sheet Model (CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glimmer_sparse

  ! This module used to be a wrapper for the umfpack and pardiso solvers.
  ! These have been removed, and now it is just a wrapper for the slap solver.

  use glimmer_global, only: dp
  use glimmer_sparse_type
  use glimmer_sparse_slap
  use glide_types

  implicit none

  type sparse_solver_options
        type(sparse_solver_options_base) :: base
        type(slap_solver_options) :: slap
  end type

  type sparse_solver_workspace
        type(slap_solver_workspace), pointer :: slap => null()
  end type

    ! The first three options use the SLAP solver and work only on one processor.
    integer, parameter :: SPARSE_SOLVER_PCG_INCH = HO_SPARSE_PCG_INCH     ! SLAP PCG with incomplete Cholesky preconditioner
    integer, parameter :: SPARSE_SOLVER_BICG = HO_SPARSE_BICG          ! SLAP biconjugate gradient
    integer, parameter :: SPARSE_SOLVER_GMRES = HO_SPARSE_GMRES         ! SLAP GMRES
    integer, parameter :: STANDALONE_PCG_STANDARD = HO_SPARSE_PCG_STANDARD     ! Native PCG, parallel-enabled, standard solver
    integer, parameter :: STANDALONE_PCG_CHRONGEAR = HO_SPARSE_PCG_CHRONGEAR    ! Native PCG, parallel-enabled, Chronopoulos-Gear solver
    integer, parameter :: STANDALONE_TRILINOS_SOLVER = HO_SPARSE_TRILINOS  ! Trilinos solver

contains



    subroutine sparse_solver_default_options(method, opt)

        integer, intent(in) :: method                 ! sparse solver: BiCG, GMRES, PCG, etc.
        type(sparse_solver_options), target :: opt    !TODO - intent inout or out?

        opt%base%method = method
        opt%base%tolerance  = 1.0d-08   !WHL - used to be 1e-11
        opt%base%maxiters = 200

        !TODO - Remove calls to not_parallel?
        !       These seem unnecessary when running SLAP solver.  Commented out for now.

        !TODO - Remove calls to slap_default_options; set appropriate options here instead.

        !Solver specific options

        if (method == SPARSE_SOLVER_BICG) then
!            call not_parallel(__FILE__,__LINE__)
            call slap_default_options(opt%slap, opt%base)
            opt%base%method = SPARSE_SOLVER_BICG
!            opt%slap%itol = 2   ! current default = 2 in slap_default_options
   
        else if (method == SPARSE_SOLVER_GMRES) then
!           call not_parallel(__FILE__,__LINE__)
            call slap_default_options(opt%slap, opt%base)
            opt%base%method = SPARSE_SOLVER_GMRES
!            opt%slap%itol = 2   ! current default = 2 in slap_default_options
             
        else if (method == SPARSE_SOLVER_PCG_INCH) then
!            call not_parallel(__FILE__, __LINE__)
            call slap_default_options(opt%slap, opt%base)
            opt%base%method = SPARSE_SOLVER_PCG_INCH
            opt%slap%itol = 1
            !WHL - itol = 2 does not work for simple test problems 

        else 
            !call glide_finalise_all(.true.)
            call write_log("Invalid sparse matrix option", GM_FATAL)

        end if

    end subroutine sparse_solver_default_options

    subroutine sparse_allocate_workspace(matrix, options, workspace, max_nonzeros_arg)

        !> Allocate solver workspace.  This needs to be done once
        !> (when the maximum number of nonzero entries is first known)
        !> This function need not be safe to call on already allocated memory
        !>
        !> Note that the max_nonzeros argument must be optional, and if
        !> it is not supplied the current number of nonzeroes must be used.
        type(sparse_matrix_type), intent(in) :: matrix
        type(sparse_solver_options) :: options
        type(sparse_solver_workspace) :: workspace
        integer, optional :: max_nonzeros_arg
        integer :: max_nonzeros
        
        if (present(max_nonzeros_arg)) then
            max_nonzeros = max_nonzeros_arg
        else
            max_nonzeros = matrix%nonzeros
        end if

        !TODO - Anything needed for standalone_pcg_solver?

        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES .or. &
            options%base%method == SPARSE_SOLVER_PCG_INCH) then
!           call not_parallel(__FILE__,__LINE__)
            allocate(workspace%slap)
            call slap_allocate_workspace(matrix, options%slap, workspace%slap, max_nonzeros)
        end if

    end subroutine sparse_allocate_workspace

    subroutine sparse_solver_preprocess(matrix, options, workspace)
        !> Performs any preprocessing needed to be performed on the slap
        !> matrix.  Workspace must have already been allocated. 
        !> This function should be safe to call more than once.
        !>
        !> It is an error to call this function on a workspace without
        !> allocated memory
        !>
        !> In general slap_allocate_workspace should perform any actions
        !> that depend on the *size* of the slap matrix, and
        !> sprase_solver_preprocess should perform any actions that depend
        !> upon the *contents* of the slap matrix.
        type(sparse_matrix_type), intent(in) :: matrix
        type(sparse_solver_options) :: options
        type(sparse_solver_workspace) :: workspace

        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES .or. &
            options%base%method == SPARSE_SOLVER_PCG_INCH) then

            call slap_solver_preprocess(matrix, options%slap, workspace%slap)

        end if

    end subroutine sparse_solver_preprocess

    function sparse_solve(matrix, rhs, solution, &
                          options, workspace,     &
                          err, niters, verbose)

        !> Solves the linear system, and reports status information.
        !> This function returns an error code that should be zero if the
        !> call succeeded and nonzero if it failed.  No additional error codes
        !> are defined.  Although this function reports back the final error
        !> and the number of iterations needed to converge, these should *not*
        !> be relied upon as not every slap linear solver may report them.

        ! Note: The matrix needs to be intent(in), not (inout).
        !      If the matrix is modified, then the residual will be computed incorrectly
        !       in the higher-level subroutine that calls sparse_solve.

        type(sparse_matrix_type), intent(in) :: matrix 
        !> Sparse matrix to solve  
        
        real(kind=dp), dimension(:), intent(in) :: rhs 
        !> Right hand side of the solution vector
        
        real(kind=dp), dimension(:), intent(inout) :: solution 
        !> Solution vector, containing an initial guess

        type(sparse_solver_options), intent(in) :: options
        !> Options such as convergence criteria
        
        type(sparse_solver_workspace), intent(inout) :: workspace
        !> Internal solver workspace
        
        real(kind=dp), intent(out) :: err
        !> Final solution error
        
        integer, intent(out) :: niters
        !> Number of iterations required to reach the solution

        logical, optional, intent(in) :: verbose
        !> If present and true, this argument may cause diagnostic information
        !> to be printed by the solver (not every solver may implement this).
        
        integer :: sparse_solve

        logical :: verbose_var

        verbose_var = .false.
        if (present(verbose)) then
            verbose_var = verbose
        end if

        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES .or. &
            options%base%method == SPARSE_SOLVER_PCG_INCH) then

            sparse_solve = slap_solve(matrix, rhs, solution, &
                                      options%slap, workspace%slap, &
                                      err, niters, verbose_var)

        end if

    end function sparse_solve


    subroutine sparse_solver_postprocess(matrix, options, workspace)
        type(sparse_matrix_type), intent(in) :: matrix
        type(sparse_solver_options) :: options
        type(sparse_solver_workspace) :: workspace

        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES .or. &
            options%base%method == SPARSE_SOLVER_PCG_INCH) then

            call slap_solver_postprocess(matrix, options%slap, workspace%slap)

        end if

    end subroutine sparse_solver_postprocess

    subroutine sparse_destroy_workspace(matrix, options, workspace)

        !> Deallocates all working memory for the slap linear solver.
        !> This need *not* be safe to call of an unallocated workspace
        !> No slap solver should call this automatically.

        type(sparse_matrix_type), intent(in) :: matrix
        type(sparse_solver_options) :: options
        type(sparse_solver_workspace) :: workspace
        
        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES .or. &
            options%base%method == SPARSE_SOLVER_PCG_INCH) then

            call slap_destroy_workspace(matrix, options%slap, workspace%slap)
            deallocate(workspace%slap)


        end if

    end subroutine sparse_destroy_workspace

    subroutine sparse_interpret_error(options, error_code, error_string)

        !> takes an error code output from slap_solve and interprets it.
        !> error_string must be an optional argument.
        !> If it is not provided, the error is printed to standard out
        !> instead of being put in the string

        type(sparse_solver_options) :: options
        integer :: error_code
        character(*), optional, intent(out) :: error_string
        character(256) :: tmp_error_string
        
        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES .or. &
            options%base%method == SPARSE_SOLVER_PCG_INCH) then

            call slap_interpret_error(error_code, tmp_error_string)

        endif

        if (present(error_string)) then
            error_string = tmp_error_string
        else
            write(*,*) tmp_error_string
        endif

    end subroutine sparse_interpret_error

    subroutine sparse_easy_solve(matrix, rhs, answer, err, iter, method_arg, &
                                 calling_file, calling_line)

        !This subroutine wraps the basic (though probably the most inefficient)
        !workflow to solve a sparse matrix using the sparse matrix solver
        !framework.  It handles errors gracefully, and reports back the
        !iterations required and the error estimate in the case of an iterative
        !solver.  At the very least it is an encapsulated example of how to
        !use the sparse solver routines, and is easy enough to drop in your
        !code if you don't care about allocating and deallocating workspace
        !every single timestep.

        type(sparse_matrix_type), intent(in) :: matrix
        real(dp), dimension(:), intent(in) :: rhs
        real(dp), dimension(:), intent(inout) :: answer
        
        real(dp), intent(out) :: err
        integer, intent(out) :: iter
        
        integer, optional, intent(in) :: method_arg         ! solver method: BiCG, GMRES, PCG, etc.

        character(*), optional :: calling_file
        integer, optional :: calling_line

        type(sparse_solver_options), target :: opt
        type(sparse_solver_workspace) :: wk

        integer :: ierr
        integer :: method

        if (present(method_arg)) then
            method = method_arg
        else
            method = SPARSE_SOLVER_BICG
        endif
 
        if (verbose_slap) then
           print*, ' '
           print*, 'In sparse_easy_solve'
           print*, 'method (0=BiCG, 1=GMRES, 2=PCG_INCH) =', method
           print*, 'matrix%order =', matrix%order
           print*, 'matrix%nonzeros =', matrix%nonzeros
           print*, 'size(rhs) =', size(rhs)
           print*, 'size(answer) =', size(answer)
           print*, 'size(row) =', size(matrix%row)
           print*, 'size(col) =', size(matrix%col)
           print*, 'size(val) =', size(matrix%val)
        endif

        call sparse_solver_default_options(method, opt)

        call sparse_allocate_workspace(matrix, opt, wk)

        call sparse_solver_preprocess(matrix, opt, wk)

        ierr = sparse_solve(matrix, rhs, answer, opt, wk, err, iter, .false.)

       if (verbose_slap) then
          print*, ' '
          print*, 'Called sparse_solve: iter, err =', iter, err
       endif
       
        call sparse_solver_postprocess(matrix, opt, wk)

        if (ierr /= 0) then
            if (present(calling_file) .and. present(calling_line)) then
                call handle_sparse_error(matrix, opt, ierr, calling_file, calling_line)
            else
                call handle_sparse_error(matrix, opt, ierr, __FILE__, __LINE__)
            end if
        end if

        call sparse_destroy_workspace(matrix, opt, wk)

    end subroutine sparse_easy_solve

    subroutine handle_sparse_error(matrix, solver_options, error, error_file, error_line, time)

        !Checks a sparse error flag and, if an error occurred, log it to
        !the GLIMMER log file.  This does not stop Glimmer, it just writes
        !to the log
        !use glide_stop
        use glimmer_log
        use glimmer_filenames
        
        integer :: error
        integer, optional :: error_line
        character(*), optional :: error_file
        real(dp), optional :: time

        type(sparse_matrix_type), intent(in) :: matrix
        type(sparse_solver_options) :: solver_options
        integer :: isym
        integer :: lunit
        integer :: i

        character(512) :: message
        character(128) :: errfname
        character(256) :: errdesc

        !If no error happened, this routine should be a nop
        if (error == 0 .OR. error == 2 .OR. error == 6) return

        !Aquire a file unit, and open the file
        lunit = get_free_unit()
        errfname = trim(process_path('sparse_dump.txt'))
        open(lunit,file=errfname)

        if (matrix%symmetric) then
            isym = 1
        else
            isym = 0
        end if

        !Output sparse matrix data to the file

        call dcpplt(matrix%order, matrix%nonzeros, matrix%row, matrix%col, matrix%val,&
                    isym, lunit)

        write(lunit,*) '***Sparse matrix structure ends.  Value listing begins'
        do i=1,matrix%nonzeros
            write(lunit,*) matrix%val(i)
        end do

        !Close unit and finish off
        close(lunit)
        
        !Grab the error message from the sparse solver
        call sparse_interpret_error(solver_options, error, errdesc)

        !construct the error message and write it to the log file
        if (present(time)) then
            write(message, *)'Sparse matrix error at time: ', time, &
                             'Error description: ', errdesc, &
                             'Data dumped to ', trim(errfname)
        else
            write(message, *)'Sparse matrix error. Error description: ', errdesc, &
                             'Data dumped to ', trim(errfname)
        end if

        write(*,*)message

        !call glide_finalise_all(.true.)

        if (present(error_file) .and. present(error_line)) then
            call write_log(trim(errdesc), GM_FATAL, error_file, error_line)
        else
            call write_log(trim(errdesc), GM_FATAL, __FILE__, __LINE__)
        end if

    end subroutine handle_sparse_error

end module glimmer_sparse
