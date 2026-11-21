!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_parallel_reduce.F90 - part of the Community Ice Sheet Model (CISM)  
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

module cism_parallel_reduce

  ! Part of the cism_parallel module:
  !  global sums, reductions (sum, max, min, maxloc, minloc) and tests for zero
  ! Other modules use these through cism_parallel.

  use glimmer_global, only : dp, sp, i8
  use glimmer_paramets, only: iulog
  use cism_reprosum_mod, only: cism_reprosum_calc, verbose_reprosum
  use cism_parallel_core

  implicit none

  interface parallel_global_sum
     module procedure parallel_global_sum_integer_2d
     module procedure parallel_global_sum_integer_3d
     module procedure parallel_global_sum_integer8_2d
     module procedure parallel_global_sum_real8_2d
     module procedure parallel_global_sum_real8_3d
  end interface

  interface parallel_global_sum_patch
     module procedure parallel_global_sum_patch_integer_2d
     module procedure parallel_global_sum_patch_real8_2d
  end interface parallel_global_sum_patch

  interface parallel_global_sum_stagger
     module procedure parallel_global_sum_stagger_integer_2d
     module procedure parallel_global_sum_stagger_real8_2d
     module procedure parallel_global_sum_stagger_real8_3d
     module procedure parallel_global_sum_stagger_real8_2d_nflds
     module procedure parallel_global_sum_stagger_real8_3d_nflds
  end interface parallel_global_sum_stagger

  interface parallel_is_zero
     module procedure parallel_is_zero_integer_1d
     module procedure parallel_is_zero_real8_1d
     module procedure parallel_is_zero_integer_2d
     module procedure parallel_is_zero_real8_2d
     module procedure parallel_is_zero_real8_3d
  end interface

  interface parallel_reduce_max
     module procedure parallel_reduce_max_integer
     module procedure parallel_reduce_max_real4
     module procedure parallel_reduce_max_real8
     module procedure parallel_reduce_max_real8_1d
  end interface

  ! This reduce interface determines the global max value and the processor on which it occurs
  interface parallel_reduce_maxloc
     module procedure parallel_reduce_maxloc_integer
     module procedure parallel_reduce_maxloc_real4
     module procedure parallel_reduce_maxloc_real8
  end interface

  interface parallel_reduce_min
     module procedure parallel_reduce_min_integer
     module procedure parallel_reduce_min_real4
     module procedure parallel_reduce_min_real8
     module procedure parallel_reduce_min_real8_1d
  end interface

  ! This reduce interface determines the global min value and the processor on which it occurs
  interface parallel_reduce_minloc
     module procedure parallel_reduce_minloc_integer
     module procedure parallel_reduce_minloc_real4
     module procedure parallel_reduce_minloc_real8
  end interface

  interface parallel_reduce_sum
     module procedure parallel_reduce_sum_integer
     module procedure parallel_reduce_sum_integer8
     module procedure parallel_reduce_sum_real4
     module procedure parallel_reduce_sum_real8
     module procedure parallel_reduce_sum_integer_nvar
     module procedure parallel_reduce_sum_real8_nvar
  end interface

contains

!=======================================================================

  ! functions belonging to the parallel_global_sum interface

  function parallel_global_sum_integer_2d(a, parallel, mask_2d)

    ! Calculates the global sum of a 2D integer field

    integer,dimension(:,:),intent(in) :: a
    type(parallel_type) :: parallel
    integer, dimension(:,:), intent(in), optional :: mask_2d

    integer :: i, j
    integer, dimension(parallel%local_ewn,parallel%local_nsn) :: mask
    integer :: local_sum
    integer :: parallel_global_sum_integer_2d

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn)

    if (present(mask_2d)) then
       mask = mask_2d
    else
       mask = 1
    endif

    local_sum = 0
    do j = nhalo+1, local_nsn-nhalo
       do i = nhalo+1, local_ewn-nhalo
          if (mask(i,j) == 1) then
             local_sum = local_sum + a(i,j)
          endif
       enddo
    enddo
    parallel_global_sum_integer_2d = parallel_reduce_sum(local_sum)

    end associate

  end function parallel_global_sum_integer_2d

!=======================================================================

  function parallel_global_sum_integer_3d(a, parallel, mask_3d)

    ! Calculates the global sum of a 3D integer field
    ! Note: The vertical dimension should be the first dimension of the input field.

    integer,dimension(:,:,:),intent(in) :: a
    type(parallel_type) :: parallel
    integer, dimension(:,:,:), intent(in), optional :: mask_3d

    integer :: i, j, k
    integer :: nz
    integer, dimension(size(a,1),parallel%local_ewn,parallel%local_nsn) :: mask
    integer :: local_sum
    integer :: parallel_global_sum_integer_3d

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn)

    nz = size(a,1)

    if (present(mask_3d)) then
       mask = mask_3d
    else
       mask = 1
    endif

    local_sum = 0
    do j = nhalo+1, local_nsn-nhalo
       do i = nhalo+1, local_ewn-nhalo
          do k = 1, nz
             if (mask(k,i,j) == 1) then
                local_sum = local_sum + a(k,i,j)
             endif
          enddo
       enddo
    enddo
    parallel_global_sum_integer_3d = parallel_reduce_sum(local_sum)

    end associate

  end function parallel_global_sum_integer_3d

!=======================================================================

  function parallel_global_sum_integer8_2d(a, parallel, mask_2d)

    ! Calculates the global sum of a 2D integer(i8) field

    integer(i8), dimension(:,:),intent(in) :: a
    type(parallel_type) :: parallel
    integer, dimension(:,:), intent(in), optional :: mask_2d

    integer :: i, j
    integer, dimension(parallel%local_ewn,parallel%local_nsn) :: mask
    integer(i8) :: local_sum
    integer(i8) :: parallel_global_sum_integer8_2d

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn)

   if (present(mask_2d)) then
       mask = mask_2d
    else
       mask = 1
    endif

    local_sum = 0
    do j = nhalo+1, local_nsn-nhalo
       do i = nhalo+1, local_ewn-nhalo
          if (mask(i,j) == 1) then
             local_sum = local_sum + a(i,j)
          endif
       enddo
    enddo
    parallel_global_sum_integer8_2d = parallel_reduce_sum(local_sum)

    end associate

  end function parallel_global_sum_integer8_2d

!=======================================================================

  function parallel_global_sum_real8_2d(a, parallel, mask_2d)

    ! Calculates the global sum of a 2D double-precision field

    real(dp), dimension(:,:), intent(in) :: a
    type(parallel_type) :: parallel
    integer, dimension(:,:), intent(in), optional :: mask_2d

    integer :: i, j
    integer, dimension(parallel%local_ewn,parallel%local_nsn) :: mask
    real(dp) :: local_sum
    real(dp) :: parallel_global_sum_real8_2d

    ! variables for computing reproductible sums
    integer :: nsummands, nflds      ! dimensions of array passed to parallel_reduce_reprosum
    integer :: count
    real(dp), dimension(:,:), allocatable :: arr
    real(dp), dimension(:), allocatable :: arr_gsum

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn)

    if (present(mask_2d)) then
       mask = mask_2d
    else
       mask = 1
    endif

    if (parallel%reprosum) then   ! compute using parallel_reduce_reprosum

       ! Allocate and fill arrays to pass to parallel_reduce_reprosum
       nsummands = own_ewn*own_nsn
       nflds = 1
       allocate(arr(nsummands,nflds))
       allocate(arr_gsum(nflds))

       count = 0
       do j = nhalo+1, local_nsn-nhalo
          do i = nhalo+1, local_ewn-nhalo
             count = count + 1
             if (mask(i,j) == 1) then
                arr(count,1) = a(i,j)
             else
                arr(count,1) = 0.0d0
             endif
          enddo
       enddo

       ! bug check
       if (count /= nsummands) then
          if (main_task) write(iulog,*) 'Error: count, nsummands =', count, nsummands
          call parallel_stop(__FILE__,__LINE__)
       endif

       ! Call parallel_reduce_reprosum
       call parallel_reduce_reprosum(arr, arr_gsum)

       parallel_global_sum_real8_2d = arr_gsum(1)

       if (verbose_reprosum .and. main_task) then
!          write(iulog,*) 'arr_gsum =', arr_gsum
       endif

       deallocate(arr)
       deallocate(arr_gsum)

    else  ! compute using parallel_reduce_sum (not reproducible)

       local_sum = 0.0d0
       do j = nhalo+1, local_nsn-nhalo
          do i = nhalo+1, local_ewn-nhalo
             if (mask(i,j) == 1) then
                local_sum = local_sum + a(i,j)
             endif
          enddo
       enddo

       ! Compute the global sum
       parallel_global_sum_real8_2d = parallel_reduce_sum(local_sum)

    endif   ! reprosum

    end associate

  end function parallel_global_sum_real8_2d

!=======================================================================

  function parallel_global_sum_real8_3d(a, parallel, mask_2d)

    ! Calculates the global sum of a 3D double-precision field
    ! Note: The vertical dimension should be the first dimension of the input field.

    real(dp), dimension(:,:,:),intent(in) :: a
    type(parallel_type) :: parallel
    integer, dimension(:,:), intent(in), optional :: mask_2d

    integer :: i, j, k
    integer :: nz
    integer, dimension(parallel%local_ewn,parallel%local_nsn) :: mask
    real(dp) :: local_sum
    real(dp) :: parallel_global_sum_real8_3d

    ! variables for computing reproductible sums
    integer :: nsummands, nflds      ! dimensions of array passed to parallel_reduce_reprosum
    integer :: count
    real(dp), dimension(:,:), allocatable :: arr
    real(dp), dimension(:), allocatable :: arr_gsum

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn)

    nz = size(a,1)

    ! Note: The mask is 2D, since typically all layers in a column are either masked in or masked out
    if (present(mask_2d)) then
       mask = mask_2d
    else
       mask = 1
    endif

    if (parallel%reprosum) then   ! compute using cism_reprosum_calc

       ! Allocate and fill arrays to pass to parallel_reduce_reprosum
       nsummands = own_ewn*own_nsn*nz
       nflds = 1
       allocate(arr(nsummands,nflds))
       allocate(arr_gsum(nflds))

       count = 0
       do j = nhalo+1, local_nsn-nhalo
          do i = nhalo+1, local_ewn-nhalo
             if (mask(i,j) == 1) then
                do k = 1, nz
                   count = count + 1
                   arr(count,1) = a(k,i,j)
                enddo
             else
                do k = 1, nz
                   count = count + 1
                   arr(count,1) = 0.0d0
                enddo
             endif
          enddo
       enddo

       ! bug check
       if (count /= nsummands) then
          if (main_task) write(iulog,*) 'Error: count, nsummands =', count, nsummands
          call parallel_stop(__FILE__,__LINE__)
       endif

       ! Call parallel_reduce_reprosum
       call parallel_reduce_reprosum(arr, arr_gsum)

       parallel_global_sum_real8_3d = arr_gsum(1)

       if (verbose_reprosum .and. main_task) then
!          write(iulog,*) 'arr_gsum =', arr_gsum
       endif

       deallocate(arr)
       deallocate(arr_gsum)

    else  ! compute using parallel_reduce_sum (not reproducible)

       local_sum = 0
       do j = nhalo+1, local_nsn-nhalo
          do i = nhalo+1, local_ewn-nhalo
             if (mask(i,j) == 1) then
                do k = 1, nz
                   local_sum = local_sum + a(k,i,j)
                enddo
             endif
          enddo
       enddo
       parallel_global_sum_real8_3d = parallel_reduce_sum(local_sum)

    endif   ! reprosum

    end associate

  end function parallel_global_sum_real8_3d

!=======================================================================
  ! subroutines belonging to the parallel_global_sum_patch interface

  function parallel_global_sum_patch_integer_2d(a, npatch, patch_id, parallel)

    ! Calculates the global sum of a 2D double-precision field over each
    !  user-defined patch of the domain.
    ! The number of patches = npatch.
    ! Each cell has an integer ID assigning it to at most one patch.
    ! If a cell has patch_id = 0, it belongs to no patches.

    integer, dimension(:,:), intent(in) :: a
    integer, intent(in) :: npatch
    integer, dimension(:,:), intent(in) :: patch_id
    type(parallel_type) :: parallel

    integer :: i, j, np
    integer, dimension(npatch) :: local_patch_sum
    integer, dimension(npatch) :: parallel_global_sum_patch_integer_2d

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn)

    local_patch_sum = 0

    do j = nhalo+1, local_nsn-nhalo
       do i = nhalo+1, local_ewn-nhalo
          np = patch_id(i,j)
          if (np > 0) then
             local_patch_sum(np) = local_patch_sum(np) + a(i,j)
          endif
       enddo
    enddo

    parallel_global_sum_patch_integer_2d = parallel_reduce_sum(local_patch_sum)

    end associate

  end function parallel_global_sum_patch_integer_2d

!=======================================================================

  function parallel_global_sum_patch_real8_2d(a, npatch, patch_id, parallel)

    ! Calculates the global sum of a 2D double-precision field over each
    !  user-defined patch of the domain.
    ! The number of patches = npatch.
    ! Each cell has an integer ID assigning it to at most one patch.
    ! If a cell has patch_id = 0, it belongs to no patches.
    !TODO - Add a reprosum option with npatch = nflds

    real(dp), dimension(:,:), intent(in) :: a
    integer, intent(in) :: npatch
    integer, dimension(:,:), intent(in) :: patch_id
    type(parallel_type) :: parallel

    integer :: i, j, np
    real(dp), dimension(npatch) :: local_patch_sum
    real(dp), dimension(npatch) :: parallel_global_sum_patch_real8_2d

    ! variables for computing reproductible sums
    integer :: nsummands, nflds      ! dimensions of array passed to parallel_reduce_reprosum
    integer :: count
    real(dp), dimension(:,:), allocatable :: arr
    real(dp), dimension(:), allocatable :: arr_gsum

    associate(  &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn)

    if (parallel%reprosum) then   ! compute using cism_reprosum_calc

       ! Allocate and fill arrays to pass to parallel_reduce_reprosum
       nsummands = own_ewn*own_nsn
       nflds = npatch
       allocate(arr(nsummands,nflds))
       allocate(arr_gsum(nflds))

       count = 0
       arr(:,:) = 0.0d0

       do j = nhalo+1, local_nsn-nhalo
          do i = nhalo+1, local_ewn-nhalo
             count = count + 1
             np = patch_id(i,j)
             if (np > 0) then
                arr(count,np) = a(i,j)
             endif
          enddo
       enddo

       ! bug check
       if (count /= nsummands) then
          if (main_task) write(iulog,*) 'Error: count, nsummands =', count, nsummands
          call parallel_stop(__FILE__,__LINE__)
       endif

       ! Call parallel_reduce_reprosum
       call parallel_reduce_reprosum(arr, arr_gsum)

       parallel_global_sum_patch_real8_2d(:) = arr_gsum(:)

       if (verbose_reprosum .and. main_task) then
!          write(iulog,*) 'arr_gsum =', arr_gsum
       endif

       deallocate(arr)
       deallocate(arr_gsum)

    else   ! compute using parallel_reduce_sum (not reproducible)

       local_patch_sum = 0.0d0

       do j = nhalo+1, local_nsn-nhalo
          do i = nhalo+1, local_ewn-nhalo
             np = patch_id(i,j)
             if (np > 0) then
                local_patch_sum(np) = local_patch_sum(np) + a(i,j)
             endif
          enddo
       enddo

       parallel_global_sum_patch_real8_2d = parallel_reduce_sum(local_patch_sum)

    endif  ! reprosum

    end associate

  end function parallel_global_sum_patch_real8_2d

!=======================================================================
  ! subroutines belonging to the parallel_global_sum_stagger interface

  function parallel_global_sum_stagger_integer_2d(a, parallel, mask_2d)

    ! Calculates the global sum of a 2D integer field on the staggered grid
    ! Similar to unstagged version, except it uses staggered_ilo/ihi/jlo/jhi

    integer,dimension(:,:),intent(in) :: a
    type(parallel_type) :: parallel
    integer, dimension(:,:), intent(in), optional :: mask_2d

    integer :: i, j
    integer, dimension(parallel%local_ewn-1,parallel%local_nsn-1) :: mask
    integer :: local_sum
    integer :: parallel_global_sum_stagger_integer_2d
    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    !TODO - associate
    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    if (present(mask_2d)) then
       mask = mask_2d
    else
       mask = 1
    endif

    local_sum = 0
    do j = staggered_jlo, staggered_jhi
       do i = staggered_ilo, staggered_ihi
          if (mask(i,j) == 1) then
             local_sum = local_sum + a(i,j)
          endif
       enddo
    enddo
    parallel_global_sum_stagger_integer_2d = parallel_reduce_sum(local_sum)

  end function parallel_global_sum_stagger_integer_2d

!=======================================================================

  function parallel_global_sum_stagger_real8_2d(arr1, parallel, arr2, mask_2d)

    ! Calculate the global sum of a 2D double-precision field on the staggered grid
    ! Similar to unstagged version, except it uses staggered_ilo/ihi/jlo/jhi

    real(dp), dimension(:,:), intent(in) :: arr1
    type(parallel_type) :: parallel
    real(dp), dimension(:,:), intent(in), optional :: arr2
    integer, dimension(:,:), intent(in), optional :: mask_2d

    integer :: i, j
    real(dp) :: local_sum
    real(dp) :: parallel_global_sum_stagger_real8_2d
    integer, dimension(parallel%local_ewn-1,parallel%local_nsn-1) :: mask

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    ! variables for computing reproductible sums
    integer :: nsummands, nflds      ! dimensions of array passed to parallel_reduce_reprosum
    integer :: count
    real(dp), dimension(:,:), allocatable :: arr
    real(dp), dimension(:), allocatable :: arr_gsum

    !TODO - associate
    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    if (present(mask_2d)) then
       mask = mask_2d
    else
       mask = 1
    endif

    if (parallel%reprosum) then   ! compute using cism_reprosum_calc

       ! Allocate and fill arrays to pass to parallel_reduce_reprosum
       nsummands = (staggered_ihi-staggered_ilo+1) * (staggered_jhi-staggered_jlo+1)
       nflds = 1
       allocate(arr(nsummands,nflds))
       allocate(arr_gsum(nflds))

       arr(:,:) = 0.0d0

       if (present(arr2)) then  ! compute global sum of arr1 + arr2

          count = 0
          do j = staggered_jlo, staggered_jhi
             do i = staggered_ilo, staggered_ihi
                count = count + 1
                if (mask(i,j) == 1) then
                   arr(count,1) = arr1(i,j) + arr2(i,j)
                endif
             enddo
          enddo

       else  ! compute global sum of arr1

          count = 0
          do j = staggered_jlo, staggered_jhi
             do i = staggered_ilo, staggered_ihi
                count = count + 1
                if (mask(i,j) == 1) then
                   arr(count,1) = arr1(i,j)
                endif
             enddo
          enddo

       endif

       ! bug check
       if (count /= nsummands) then
          if (main_task) write(iulog,*) 'Error: count, nsummands =', count, nsummands
          call parallel_stop(__FILE__,__LINE__)
       endif

       ! Call parallel_reduce_reprosum
       call parallel_reduce_reprosum(arr, arr_gsum)

       parallel_global_sum_stagger_real8_2d = arr_gsum(1)

       if (verbose_reprosum .and. main_task) then
!          write(iulog,*) 'arr_gsum =', arr_gsum
       endif

       deallocate(arr)
       deallocate(arr_gsum)

    else  ! compute using parallel_reduce_sum (not reproducible)

       local_sum = 0.0d0

       if (present(arr2)) then  ! compute global sum of arr1 + arr2

          do j = staggered_jlo, staggered_jhi
             do i = staggered_ilo, staggered_ihi
                if (mask(i,j) == 1) then
                   local_sum = local_sum + arr1(i,j) + arr2(i,j)
                endif
             enddo
          enddo

       else  ! compute global sum of arr1

          do j = staggered_jlo, staggered_jhi
             do i = staggered_ilo, staggered_ihi
                if (mask(i,j) == 1) then
                   local_sum = local_sum + arr1(i,j)
                endif
             enddo
          enddo

       endif

       parallel_global_sum_stagger_real8_2d = parallel_reduce_sum(local_sum)

    endif   ! reprosum

  end function parallel_global_sum_stagger_real8_2d

!=======================================================================

  !TODO - Add mask_2d here and below
  function parallel_global_sum_stagger_real8_3d(arr1, parallel, arr2, mask_2d)

    ! Calculate the global sum of a 3D double-precision field on the staggered grid
    ! Assumes k is the first index, followed by i and j

    real(dp), dimension(:,:,:), intent(in) :: arr1
    type(parallel_type) :: parallel
    real(dp), dimension(:,:,:), intent(in), optional :: arr2
    integer, dimension(:,:), intent(in), optional :: mask_2d

    integer :: i, j, k, nz
    real(dp) :: local_sum
    real(dp) :: parallel_global_sum_stagger_real8_3d
    integer, dimension(parallel%local_ewn-1,parallel%local_nsn-1) :: mask

    ! variables for computing reproductible sums
    integer :: nsummands, nflds      ! dimensions of array passed to parallel_reduce_reprosum
    integer :: count
    real(dp), dimension(:,:), allocatable :: arr
    real(dp), dimension(:), allocatable :: arr_gsum

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    !TODO - Associate these variables (and not the ones above)
    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    if (present(mask_2d)) then
       mask = mask_2d
    else
       mask = 1
    endif

    nz = size(arr1,1)

    if (parallel%reprosum) then   ! compute using cism_reprosum_calc

       ! Allocate and fill arrays to pass to parallel_reduce_reprosum
       nsummands = (staggered_ihi-staggered_ilo+1) * (staggered_jhi-staggered_jlo+1) * nz
       nflds = 1
       allocate(arr(nsummands,nflds))
       allocate(arr_gsum(nflds))

       arr(:,:) = 0.0d0

       if (present(arr2)) then  ! compute global sum of arr1 + arr2

          count = 0
          do j = staggered_jlo, staggered_jhi
             do i = staggered_ilo, staggered_ihi
                if (mask(i,j) == 1) then
                   do k = 1, nz
                      count = count + 1
                      arr(count,1) = arr1(k,i,j) + arr2(k,i,j)
                   enddo
                endif
             enddo
          enddo

       else  ! compute global sum of arr1

          count = 0
          do j = staggered_jlo, staggered_jhi
             do i = staggered_ilo, staggered_ihi
                if (mask(i,j) == 1) then
                   do k = 1, nz
                      count = count + 1
                      arr(count,1) = arr1(k,i,j)
                   enddo
                endif
             enddo
          enddo

       endif

       ! bug check
       if (count /= nsummands) then
          if (main_task) write(iulog,*) 'Error: count, nsummands =', count, nsummands
          call parallel_stop(__FILE__,__LINE__)
       endif

       ! Call parallel_reduce_reprosum
       call parallel_reduce_reprosum(arr, arr_gsum)

       parallel_global_sum_stagger_real8_3d = arr_gsum(1)

       if (verbose_reprosum .and. main_task) then
!          write(iulog,*) 'arr_gsum =', arr_gsum
       endif

       deallocate(arr)
       deallocate(arr_gsum)

    else  ! compute using parallel_reduce_sum (not reproducible)

       local_sum = 0.0d0

       if (present(arr2)) then  ! compute global sum of arr1 + arr2

          do j = staggered_jlo, staggered_jhi
             do i = staggered_ilo, staggered_ihi
                if (mask(i,j) == 1) then
                   do k = 1, nz
                      local_sum = local_sum + arr1(k,i,j) + arr2(k,i,j)
                   enddo
                endif
             enddo
          enddo

       else  ! compute global sum of arr1

          do j = staggered_jlo, staggered_jhi
             do i = staggered_ilo, staggered_ihi
                if (mask(i,j) == 1) then
                   do k = 1, nz
                      local_sum = local_sum + arr1(k,i,j)
                   enddo
                endif
             enddo
          enddo

       endif

       parallel_global_sum_stagger_real8_3d = parallel_reduce_sum(local_sum)

    endif   ! reprosum

  end function parallel_global_sum_stagger_real8_3d

!=======================================================================

  function parallel_global_sum_stagger_real8_2d_nflds(arr1, nflds, parallel, arr2, mask_2d)

    ! Sum one or two local arrays on the staggered grid, then take the global sum.
    ! The final index is equal to the number of independent fields to be summed.
    !TODO - Don't have to pass in nflds, since it equals size(a,3)?

    real(dp), dimension(:,:,:), intent(in) :: arr1
    integer, intent(in) :: nflds
    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication
    real(dp), dimension(:,:,:), intent(in), optional :: arr2
    integer, dimension(:,:), intent(in), optional :: mask_2d

    real(dp), dimension(size(arr1,3)) :: parallel_global_sum_stagger_real8_2d_nflds
    integer, dimension(parallel%local_ewn-1,parallel%local_nsn-1) :: mask

    integer :: i, j, n

    real(dp), dimension(size(arr1,3)) :: local_sum

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    ! variables for computing reproductible sums
    integer :: nsummands               ! dimensions of array passed to parallel_reduce_reprosum
    integer :: count
    real(dp), dimension(:,:), allocatable :: arr
    real(dp), dimension(:), allocatable :: arr_gsum

    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    if (present(mask_2d)) then
       mask = mask_2d
    else
       mask = 1
    endif

    if (parallel%reprosum) then   ! compute using cism_reprosum_calc

       ! Allocate and fill arrays to pass to parallel_reduce_reprosum
       nsummands = (staggered_ihi-staggered_ilo+1) * (staggered_jhi-staggered_jlo+1)
       allocate(arr(nsummands,nflds))
       allocate(arr_gsum(nflds))

       arr(:,:) = 0.0d0

       do n = 1, nflds

          if (present(arr2)) then  ! compute global sum of arr1 + arr2
             count = 0
             do j = staggered_jlo, staggered_jhi
                do i = staggered_ilo, staggered_ihi
                   if (mask(i,j) == 1) then
                      count = count + 1
                      arr(count,n) = arr1(i,j,n) + arr2(i,j,n)
                   endif
                enddo
             enddo
          else  ! compute global sum of arr1
             count = 0
             do j = staggered_jlo, staggered_jhi
                do i = staggered_ilo, staggered_ihi
                   if (mask(i,j) == 1) then
                      count = count + 1
                      arr(count,n) = arr1(i,j,n)
                   endif
                enddo
             enddo
          endif

       enddo   ! nflds

       ! bug check
       if (count /= nsummands) then
          if (main_task) write(iulog,*) 'Error: count, nsummands =', count, nsummands
          call parallel_stop(__FILE__,__LINE__)
       endif

       ! Call parallel_reduce_reprosum
       call parallel_reduce_reprosum(arr, arr_gsum)

       parallel_global_sum_stagger_real8_2d_nflds = arr_gsum(:)

       if (verbose_reprosum .and. main_task) then
!          write(iulog,*) 'arr_gsum =', arr_gsum
       endif

       deallocate(arr)
       deallocate(arr_gsum)

    else  ! compute using parallel_reduce_sum (not reproducible)

       local_sum(:) = 0.d0

       do n = 1, nflds

          ! sum over locally owned velocity points

          if (present(arr2)) then
             do j = staggered_jlo, staggered_jhi
                do i = staggered_ilo, staggered_ihi
                   if (mask(i,j) == 1) then
                      local_sum(n) = local_sum(n) + arr1(i,j,n) + arr2(i,j,n)
                   endif
                enddo
             enddo
          else
             do j = staggered_jlo, staggered_jhi
                do i = staggered_ilo, staggered_ihi
                   if (mask(i,j) == 1) then
                      local_sum(n) = local_sum(n) + arr1(i,j,n)
                   endif
                enddo
             enddo
          endif

       enddo   ! nflds

       parallel_global_sum_stagger_real8_2d_nflds = parallel_reduce_sum(local_sum(:))

    endif   ! reprosum

  end function parallel_global_sum_stagger_real8_2d_nflds

!=======================================================================

  function parallel_global_sum_stagger_real8_3d_nflds(arr1, nflds, parallel, arr2, mask_2d)

    ! Sum one or two local arrays on the staggered grid, then take the global sum.
    ! Assumes k is the first index, followed by i and j.
    ! The final index is equal to the number of independent fields to be summed.

    real(dp), dimension(:,:,:,:), intent(in) :: arr1
    integer, intent(in) :: nflds   ! size of final index; number of global sums to be computed
    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication
    real(dp), dimension(:,:,:,:), intent(in), optional :: arr2
    integer, dimension(:,:), intent(in), optional :: mask_2d

    real(dp), dimension(size(arr1,4)) :: parallel_global_sum_stagger_real8_3d_nflds
    integer, dimension(parallel%local_ewn-1,parallel%local_nsn-1) :: mask

    integer :: i, j, k, n, nz

    real(dp), dimension(size(arr1,4)) :: local_sum

    integer :: &
         staggered_ilo, staggered_ihi, &  ! bounds of locally owned vertices on staggered grid
         staggered_jlo, staggered_jhi

    ! variables for computing reproductible sums
    integer :: nsummands               ! dimensions of array passed to parallel_reduce_reprosum
    integer :: count
    real(dp), dimension(:,:), allocatable :: arr
    real(dp), dimension(:), allocatable :: arr_gsum

    !TODO - Associate these variables
    staggered_ilo = parallel%staggered_ilo
    staggered_ihi = parallel%staggered_ihi
    staggered_jlo = parallel%staggered_jlo
    staggered_jhi = parallel%staggered_jhi

    if (present(mask_2d)) then
       mask = mask_2d
    else
       mask = 1
    endif

    nz = size(arr1,1)

    if (parallel%reprosum) then   ! compute using cism_reprosum_calc

       ! Allocate and fill arrays to pass to parallel_reduce_reprosum
       nsummands = (staggered_ihi-staggered_ilo+1) * (staggered_jhi-staggered_jlo+1) * nz
       allocate(arr(nsummands,nflds))
       allocate(arr_gsum(nflds))

       arr(:,:) = 0.0d0

       do n = 1, nflds

          if (present(arr2)) then  ! compute global sum of arr1 + arr2
             count = 0
             do j = staggered_jlo, staggered_jhi
                do i = staggered_ilo, staggered_ihi
                   if (mask(i,j) == 1) then
                      do k = 1, nz
                         count = count + 1
                         arr(count,n) = arr1(k,i,j,n) + arr2(k,i,j,n)
                      enddo
                   endif
                enddo
             enddo
          else  ! compute global sum of arr1
             count = 0
             do j = staggered_jlo, staggered_jhi
                do i = staggered_ilo, staggered_ihi
                   if (mask(i,j) == 1) then
                      do k = 1, nz
                         count = count + 1
                         arr(count,n) = arr1(k,i,j,n)
                      enddo
                   endif
                enddo
             enddo
          endif

       enddo   ! nflds

       ! bug check
       if (count /= nsummands) then
          if (main_task) write(iulog,*) 'Error: count, nsummands =', count, nsummands
          call parallel_stop(__FILE__,__LINE__)
       endif

       ! Call parallel_reduce_reprosum
       call parallel_reduce_reprosum(arr, arr_gsum)

       parallel_global_sum_stagger_real8_3d_nflds = arr_gsum(:)

       if (verbose_reprosum .and. main_task) then
!          write(iulog,*) 'arr_gsum =', arr_gsum
       endif

       deallocate(arr)
       deallocate(arr_gsum)

    else  ! compute using parallel_reduce_sum (not reproducible)

       local_sum(:) = 0.d0

       do n = 1, nflds

          ! sum over locally owned velocity points

          if (present(arr2)) then
             do j = staggered_jlo, staggered_jhi
                do i = staggered_ilo, staggered_ihi
                   if (mask(i,j) == 1) then
                      do k = 1, nz
                         local_sum(n) = local_sum(n) + arr1(k,i,j,n) + arr2(k,i,j,n)
                      enddo
                   endif
                enddo
             enddo
          else
             do j = staggered_jlo, staggered_jhi
                do i = staggered_ilo, staggered_ihi
                   if (mask(i,j) == 1) then
                      do k = 1, nz
                         local_sum(n) = local_sum(n) + arr1(k,i,j,n)
                      enddo
                   endif
                enddo
             enddo
          endif

       enddo   ! nflds

       parallel_global_sum_stagger_real8_3d_nflds = parallel_reduce_sum(local_sum(:))

    endif   ! reprosum

  end function parallel_global_sum_stagger_real8_3d_nflds

!=======================================================================
  ! functions belonging to the parallel_is_zero interface

  function parallel_is_zero_integer_1d(a)

    ! returns .true. if the field has all zero values, else returns .false.

    integer, dimension(:), intent(in) :: a
    logical :: parallel_is_zero_integer_1d

    integer :: maxval_a

    maxval_a = maxval(abs(a))
    maxval_a = parallel_reduce_max(maxval_a)
    if (maxval_a > 0) then
       parallel_is_zero_integer_1d = .false.
    else
       parallel_is_zero_integer_1d = .true.
    endif

  end function parallel_is_zero_integer_1d

!=======================================================================

  function parallel_is_zero_real8_1d(a)

    ! returns .true. if the field has all zero values, else returns .false.

    real(dp), dimension(:), intent(in) :: a
    logical :: parallel_is_zero_real8_1d

    real(dp) :: maxval_a

    maxval_a = maxval(abs(a))
    maxval_a = parallel_reduce_max(maxval_a)
    if (maxval_a > 0.0d0) then
       parallel_is_zero_real8_1d = .false.
    else
       parallel_is_zero_real8_1d = .true.
    endif

  end function parallel_is_zero_real8_1d

!=======================================================================

  function parallel_is_zero_integer_2d(a)

    ! returns .true. if the field has all zero values, else returns .false.

    integer, dimension(:,:), intent(in) :: a
    logical :: parallel_is_zero_integer_2d

    integer :: maxval_a

    maxval_a = maxval(abs(a))
    maxval_a = parallel_reduce_max(maxval_a)
    if (maxval_a > 0) then
       parallel_is_zero_integer_2d = .false.
    else
       parallel_is_zero_integer_2d = .true.
    endif

  end function parallel_is_zero_integer_2d

!=======================================================================

  function parallel_is_zero_real8_2d(a)

    ! returns .true. if the field has all zero values, else returns .false.

    real(dp), dimension(:,:), intent(in) :: a
    logical :: parallel_is_zero_real8_2d

    real(dp) :: maxval_a

    maxval_a = maxval(abs(a))
    maxval_a = parallel_reduce_max(maxval_a)
    if (maxval_a > 0.0d0) then
       parallel_is_zero_real8_2d = .false.
    else
       parallel_is_zero_real8_2d = .true.
    endif

  end function parallel_is_zero_real8_2d

!=======================================================================

  function parallel_is_zero_real8_3d(a)

    ! returns .true. if the field has all zero values, else returns .false.

    real(dp), dimension(:,:,:), intent(in) :: a
    logical :: parallel_is_zero_real8_3d

    real(dp) :: maxval_a

    maxval_a = maxval(abs(a))
    maxval_a = parallel_reduce_max(maxval_a)
    if (maxval_a > 0.0d0) then
       parallel_is_zero_real8_3d = .false.
    else
       parallel_is_zero_real8_3d = .true.
    endif

  end function parallel_is_zero_real8_3d

!=======================================================================

  ! functions for parallel reduction of logical variables
  ! * parallel_reduce_log_or returns 'true' iff x = 'true' on at least one processor
  ! * parallel_reduce_log_and returns 'true' iff x = 'true' on all processors

  function parallel_reduce_log_or(x)

    use mpi_mod
    implicit none
    logical :: x

    integer :: ierror
    logical :: recvbuf,sendbuf, parallel_reduce_log_or

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_logical,mpi_lor,comm,ierror)
    parallel_reduce_log_or = recvbuf

  end function parallel_reduce_log_or

  function parallel_reduce_log_and(x)

    use mpi_mod
    implicit none
    logical :: x

    integer :: ierror
    logical :: recvbuf,sendbuf, parallel_reduce_log_and

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_logical,mpi_land,comm,ierror)
    parallel_reduce_log_and = recvbuf

  end function parallel_reduce_log_and

!=======================================================================

  ! functions belonging to the parallel_reduce_sum interface

  function parallel_reduce_sum_integer(x)

    use mpi_mod
    implicit none
    integer :: x

    integer :: ierror
    integer :: recvbuf,sendbuf, parallel_reduce_sum_integer

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_integer,mpi_sum,comm,ierror)
    parallel_reduce_sum_integer = recvbuf

  end function parallel_reduce_sum_integer


  function parallel_reduce_sum_integer8(x)

    use mpi_mod
    implicit none
    integer(i8) :: x

    integer :: ierror
    integer(i8) :: recvbuf,sendbuf, parallel_reduce_sum_integer8

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_integer8,mpi_sum,comm,ierror)
    parallel_reduce_sum_integer8 = recvbuf

  end function parallel_reduce_sum_integer8


  function parallel_reduce_sum_real4(x)

    use mpi_mod
    implicit none
    real(sp) :: x

    integer :: ierror
    real(sp) :: recvbuf,sendbuf, parallel_reduce_sum_real4

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real4,mpi_sum,comm,ierror)
    parallel_reduce_sum_real4 = recvbuf

  end function parallel_reduce_sum_real4


  function parallel_reduce_sum_real8(x)

    use mpi_mod

    implicit none

    real(dp) :: x

    integer :: ierror
    real(dp) :: recvbuf,sendbuf, parallel_reduce_sum_real8

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real8,mpi_sum,comm,ierror)
    parallel_reduce_sum_real8 = recvbuf

  end function parallel_reduce_sum_real8


  function parallel_reduce_sum_integer_nvar(x)

    use mpi_mod
    implicit none
    integer :: x(:)

    integer :: ierror, nvar
    integer, dimension(size(x)) :: recvbuf,sendbuf, parallel_reduce_sum_integer_nvar

    ! begin
    nvar = size(x)
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,nvar,mpi_integer,mpi_sum,comm,ierror)
    parallel_reduce_sum_integer_nvar = recvbuf

  end function parallel_reduce_sum_integer_nvar


  function parallel_reduce_sum_real8_nvar(x)

    use mpi_mod
    implicit none
    real(dp) :: x(:)

    integer :: ierror, nvar
    real(dp), dimension(size(x)) :: recvbuf,sendbuf, parallel_reduce_sum_real8_nvar

    ! begin
    nvar = size(x)
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,nvar,mpi_real8,mpi_sum,comm,ierror)
    parallel_reduce_sum_real8_nvar = recvbuf

  end function parallel_reduce_sum_real8_nvar

!=======================================================================

  subroutine parallel_reduce_reprosum(arr, arr_gsum)

    ! Compute a reproducible global sum for a floating-point variable or array.
    ! Can be called from parallel_global_sum, parallel_global_sum_patch, or
    ! parallel_global_sum_stagger.

    implicit none

    real(dp), dimension(:,:), intent(in) :: arr
    real(dp), dimension(:), intent(out) :: arr_gsum

    ! Notes on subroutine cism_reprosum_calc:
    ! The first five arguments are required: arr(dsummands,nflds), arr_gsum(dsummands),
    !  dsummands, dflds and nsummands. Typically, nsummands = dsummands = number of local values.
    ! We use the default fixed-precision algorithm (instead of ddpdd, which apparently is less robust).
    ! We do not allow Inf or NaN values in the input array.
    ! Typically, the algorithm calls mpi_allreduce twice. By passing in both arr_gbl_max and arr_max_levels,
    !  it may be possible to call mpi_allreduce just once, improving performance.
    !  If we don't pass these arguments, then arr_gbl_max and arr_max_levels are computed internally.
    !  By passing arr_glb_max_out and arr_max_levels_out, we can see the calculated values.
    ! By passing rel_diff, we can verify that the computed reproducible sum is close
    !  to the (nonreproducible) floating-point value.
    ! Note: In cism_reprosum_mod, the default communicator ID is MPI_COMM_WORLD.
    !       We need to pass in the CISM communicator ID ('comm') in case comm /= MPI_COMM_WORLD.
    !       CESM coupled runs fail if comm is not passed in.
    ! See comments in cism_reprosum_calc for more info

    ! Required arguments

    integer :: dsummands, nflds                      ! dimensions of arr
    integer :: nsummands                             ! number of processors

    ! Optional arguments

    real(dp), dimension(:), allocatable :: &
         arr_gbl_max,         & ! upper bound on max(abs(arr))
         arr_gbl_max_out        ! calculated upper bound on max(abs(arr))

    real(dp), dimension(:,:), allocatable :: &
         rel_diff               ! relative and absolute differences between fixed and floating point sums

    integer, dimension(:), allocatable :: &
         arr_max_levels,      & ! maximum number of levels of integer expansion to use
         arr_max_levels_out     ! output of number of levels of integer expansion to use

    integer :: &
         gbl_max_nsummands,   & ! maximum of nsummand over all processes
         gbl_max_nsummands_out  ! calculated maximum nsummands over all processes

    integer, dimension(6) :: &
         repro_sum_stats(6)     ! increment running totals for
                                !  (1) one-reduction repro_sum
                                !  (2) two-reduction repro_sum
                                !  (3) both types in one call
                                !  (4) nonrepro_sum
                                !  (5) global max nsummands reduction
                                !  (6) global lor 3*nflds reduction

    logical :: &
         ddpdd_sum,           & ! use ddpdd algorithm instead of fixed-precision algorithm
         repro_sum_validate     ! flag enabling/disabling testing that gmax and max_levels
                                ! are accurate/sufficient. Default is enabled.

    ! Set parameters and allocate arrays
    dsummands = size(arr,1)
    nflds = size(arr,2)
    nsummands = dsummands

    allocate (arr_gbl_max(nflds))
    allocate (arr_gbl_max_out(nflds))
    allocate (arr_max_levels(nflds))
    allocate (arr_max_levels_out(nflds))
    allocate (rel_diff(2,nflds))

    ddpdd_sum = .false.
    repro_sum_validate = .true.

    ! The following subroutine is adapted from shr_reprosum_calc in CESM shared code.

    call cism_reprosum_calc(&
         arr,        arr_gsum,                          &
         nsummands,  dsummands,  nflds,                 &
         ddpdd_sum = ddpdd_sum,                         &
!            arr_gbl_max = arr_gbl_max,                     &
         arr_gbl_max_out = arr_gbl_max_out,             &
!            arr_max_levels = arr_max_levels,               &
         arr_max_levels_out = arr_max_levels_out,       &
!            gbl_max_nsummands = gbl_max_nsummands,         &
         gbl_max_nsummands_out = gbl_max_nsummands_out, &
         repro_sum_validate = repro_sum_validate,       &
         repro_sum_stats = repro_sum_stats,             &
         rel_diff = rel_diff,                           &
         commid = comm)

    if (verbose_reprosum .and. main_task) then
!       write(iulog,*) 'arr_gbl_max_out =', arr_gbl_max_out
!       write(iulog,*) 'arr_max_levels_out =', arr_max_levels_out
!       write(iulog,*) 'gbl_max_nsummands_out =', gbl_max_nsummands_out
!       write(iulog,*) 'rel diff =', rel_diff(1,:)
!       write(iulog,*) 'abs diff =', rel_diff(2,:)
!       write(iulog,*) 'stats =', repro_sum_stats(:)
    endif

    deallocate(arr_gbl_max, arr_gbl_max_out)
    deallocate(arr_max_levels, arr_max_levels_out)
    deallocate(rel_diff)

  end subroutine parallel_reduce_reprosum

!=======================================================================

  ! functions belonging to the parallel_reduce_max interface

  function parallel_reduce_max_integer(x)

    use mpi_mod
    implicit none
    integer :: x

    integer :: ierror
    integer :: recvbuf,sendbuf, parallel_reduce_max_integer

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_integer,mpi_max,comm,ierror)
    parallel_reduce_max_integer = recvbuf

  end function parallel_reduce_max_integer


  function parallel_reduce_max_real4(x)

    use mpi_mod
    implicit none
    real(sp) :: x

    integer :: ierror
    real(sp) :: recvbuf,sendbuf, parallel_reduce_max_real4

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real4,mpi_max,comm,ierror)
    parallel_reduce_max_real4 = recvbuf

  end function parallel_reduce_max_real4


  function parallel_reduce_max_real8(x)

    use mpi_mod
    implicit none
    real(dp) :: x

    integer :: ierror
    real(dp) :: recvbuf,sendbuf, parallel_reduce_max_real8

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real8,mpi_max,comm,ierror)
    parallel_reduce_max_real8 = recvbuf

  end function parallel_reduce_max_real8

  function parallel_reduce_max_real8_1d(x)

    use mpi_mod
    implicit none
    real(dp), dimension(:) :: x

    integer :: ierror
    real(dp), dimension(size(x)) :: recvbuf,sendbuf, parallel_reduce_max_real8_1d

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,size(x),mpi_real8,mpi_max,comm,ierror)
    parallel_reduce_max_real8_1d = recvbuf

  end function parallel_reduce_max_real8_1d

!=======================================================================

  ! functions belonging to the parallel_reduce_maxloc interface

  subroutine parallel_reduce_maxloc_integer(xin, xout, xprocout)

    use mpi_mod
    implicit none
    integer, intent(in) :: xin         ! variable to reduce
    integer, intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    integer :: ierror
    integer, dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2INTEGER,mpi_maxloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1)

  end subroutine parallel_reduce_maxloc_integer


  subroutine parallel_reduce_maxloc_real4(xin, xout, xprocout)

    use mpi_mod
    implicit none
    real(sp), intent(in) :: xin         ! variable to reduce
    real(sp), intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout    ! processor on which reduced value occurs

    integer :: ierror
    real(sp), dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x (coerced to a real)
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2REAL,mpi_maxloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1) ! coerced back to integer

  end subroutine parallel_reduce_maxloc_real4


  subroutine parallel_reduce_maxloc_real8(xin, xout, xprocout)

    use mpi_mod
    implicit none
    real(dp), intent(in) :: xin         ! variable to reduce
    real(dp), intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout    ! processor on which reduced value occurs

    integer :: ierror
    real(dp), dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x (coerced to a real)
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2DOUBLE_PRECISION,mpi_maxloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1) ! coerced back to integer

  end subroutine parallel_reduce_maxloc_real8


!=======================================================================

  ! subroutines belonging to the parallel_reduce_min interface

  function parallel_reduce_min_integer(x)

    use mpi_mod
    implicit none
    integer :: x

    integer :: ierror
    integer :: recvbuf,sendbuf, parallel_reduce_min_integer

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_integer,mpi_min,comm,ierror)
    parallel_reduce_min_integer = recvbuf

  end function parallel_reduce_min_integer


  function parallel_reduce_min_real4(x)

    use mpi_mod
    implicit none
    real(sp) :: x

    integer :: ierror
    real(sp) :: recvbuf,sendbuf, parallel_reduce_min_real4

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real4,mpi_min,comm,ierror)
    parallel_reduce_min_real4 = recvbuf

  end function parallel_reduce_min_real4


  function parallel_reduce_min_real8(x)

    use mpi_mod
    implicit none
    real(dp) :: x

    integer :: ierror
    real(dp) :: recvbuf,sendbuf, parallel_reduce_min_real8

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,1,mpi_real8,mpi_min,comm,ierror)
    parallel_reduce_min_real8 = recvbuf

  end function parallel_reduce_min_real8


  function parallel_reduce_min_real8_1d(x)

    use mpi_mod
    implicit none
    real(dp), dimension(:) :: x

    integer :: ierror
    real(dp), dimension(size(x)) :: recvbuf,sendbuf, parallel_reduce_min_real8_1d

    ! begin
    sendbuf = x
    call mpi_allreduce(sendbuf,recvbuf,size(x),mpi_real8,mpi_min,comm,ierror)
    parallel_reduce_min_real8_1d = recvbuf

  end function parallel_reduce_min_real8_1d

!=======================================================================

  ! subroutines belonging to the parallel_reduce_minloc interface

  subroutine parallel_reduce_minloc_integer(xin, xout, xprocout)

    use mpi_mod
    implicit none
    integer, intent(in) :: xin         ! variable to reduce
    integer, intent(out) :: xout       ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    integer :: ierror
    integer, dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2INTEGER,mpi_minloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1)

  end subroutine parallel_reduce_minloc_integer


  subroutine parallel_reduce_minloc_real4(xin, xout, xprocout)

    use mpi_mod
    implicit none
    real(sp), intent(in) :: xin        ! variable to reduce
    real(sp), intent(out) :: xout      ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    integer :: ierror
    real(sp), dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x (coerced to a real)
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2REAL,mpi_minloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1) ! coerced back to integer

  end subroutine parallel_reduce_minloc_real4


  subroutine parallel_reduce_minloc_real8(xin, xout, xprocout)

    use mpi_mod
    implicit none
    real(dp), intent(in) :: xin        ! variable to reduce
    real(dp), intent(out) :: xout      ! value resulting from the reduction
    integer, intent(out) :: xprocout   ! processor on which reduced value occurs

    integer :: ierror
    real(dp), dimension(2,1) :: recvbuf, sendbuf

    ! begin
    sendbuf(1,1) = xin
    sendbuf(2,1) = this_rank  ! This is the processor number associated with the value x (coerced to a real)
    call mpi_allreduce(sendbuf,recvbuf,1,MPI_2DOUBLE_PRECISION,mpi_minloc,comm,ierror)
    xout = recvbuf(1,1)
    xprocout = recvbuf(2,1) ! coerced back to integer

  end subroutine parallel_reduce_minloc_real8

!=======================================================================

  subroutine parallel_show_minmax(label,values)

    use mpi_mod
    implicit none
    character(*) :: label
    real(dp),dimension(:,:,:) :: values
    
    integer :: ierror
    real(dp) :: allmin,allmax,mymin,mymax

    ! begin
    mymin = minval(values(:,1+lhalo:size(values,2)-uhalo,&
         1+lhalo:size(values,3)-uhalo))
    mymax = maxval(values(:,1+lhalo:size(values,2)-uhalo,&
         1+lhalo:size(values,3)-uhalo))
    call mpi_reduce(mymin,allmin,1,mpi_real8,mpi_min,main_rank,comm,ierror)
    call mpi_reduce(mymax,allmax,1,mpi_real8,mpi_max,main_rank,comm,ierror)
    if (main_task) print *,label,allmin,allmax

  end subroutine parallel_show_minmax

!=======================================================================

end module cism_parallel_reduce

!=======================================================================
