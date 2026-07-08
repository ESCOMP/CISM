
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_utils.F90 - part of the Community Ice Sheet Model (CISM)  
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

!> Module containing utility code for Glimmer and CISM
!TODO - Move check_conformal and fix_bcs2d to Glint?  Used by glint_interp only.

module glimmer_utils

  use glimmer_global, only: dp
  use glimmer_paramets, only: iulog

  implicit none

  interface array_bcs
     module procedure array_bcs1d,array_bcs2d
  end interface

  interface check_conformal
     module procedure check_conformal_2d_real
  end interface

  interface point_diag
     module procedure point_diag_integer_2d
     module procedure point_diag_logical_2d
     module procedure point_diag_real8_2d
  end interface point_diag

contains

  !> Returns the value of a 1D array location,checking first for the boundaries.
  !!
  !! the location is wrapped around the array boundaries until it falls within the array
  !! \author The value of the location in question.
  real(dp) function array_bcs1d(array,i)

    ! Arguments

    real(dp),dimension(:),intent(in) :: array !< The array to be indexed.
    integer,intent(in)               :: i     !< The location to be extracted.

    ! Internal variables

    integer :: n,ii

    n=size(array)
    ii=i

    if ((i<=n).and.(i>=1)) then
      array_bcs1d=array(i)
    endif

    do while (ii>n)
      ii=ii-n
    enddo

    do while (ii<1)
      ii=ii+n
    enddo

    array_bcs1d=array(ii)

  end function array_bcs1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Returns the value of a 1D array location,checking first for the boundaries.
  !!
  !! the location is wrapped around the array boundaries until it falls within the array
  !! as array_bcs1d but for polar boundary conditions
  !! \author The value of the location in question.
  real(dp) function array_bcs_lats(array,i)


    ! Arguments

    real(dp),dimension(:),intent(in) :: array !< The array to be indexed.
    integer,intent(in) :: i !< The location to be extracted.

    ! Internal variables

    integer :: n,ii

    n=size(array)
    ii=i

    if ((i<=n).and.(i>=1)) then
      array_bcs_lats=array(i)
      return
    endif

    if (ii>n) then
      ii=2*n-ii
      array_bcs_lats=-180.d0+array(ii)
    endif

    if (ii<1) then
      ii=1-ii
      array_bcs_lats=180.d0-array(ii)
    endif

  end function array_bcs_lats

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Returns the value of an array 
  !! location, checking first for the boundaries. 
  !! Over-the-pole boundary conditions are implemented here.
  !! \return The value of the location specified.
  real(dp) function array_bcs2d(array,i,j)

    ! Arguments

    real(dp),dimension(:,:),intent(in) :: array !< Array to be indexed
    integer,intent(in) :: i !< The location to be extracted    
    integer,intent(in) :: j !< The location to be extracted

    ! Internal variables

    integer :: nx,ny,ii,jj

    nx=size(array,1) ; ny=size(array,2)

    if ((i>=1).and.(i<=nx).and.(j>=1).and.(j<=ny)) then
      array_bcs2d=array(i,j)
      return
    endif

    ii=i ; jj=j

    if (jj>ny) then
      jj=2*ny-jj
      ii=ii+nx/2
    endif

    if (jj<1) then
      jj=1-jj
      ii=ii+nx/2
    endif

    do while (ii>nx) 
      ii=ii-nx
    enddo

    do while (ii<1)
      ii=ii+nx
    enddo  

    array_bcs2d=array(ii,jj)

  end function array_bcs2d

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine fix_bcs2d(i,j,nx,ny)

  !> Adjusts array location indices
  !! so that they fall within the domain.

    integer,intent(inout) :: i !< The location of interest
    integer,intent(inout) :: j !< The location of interest
    integer,intent(in) :: nx  !< The size of the domain (number of points in each direction)
    integer,intent(in) :: ny  !< The size of the domain (number of points in each direction)

    if ((i>=1).and.(i<=nx).and.(j>=1).and.(j<=ny)) return

    if (j>ny) then
      j=2*ny-j
      i=i+nx/2
    endif

    if (j<1) then
      j=1-j
      i=i+nx/2
    endif

    do while (i>nx) 
      i=i-nx
    enddo

    do while (i<1)
      i=i+nx
    enddo  

  end subroutine fix_bcs2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check_conformal_2d_real(array1,array2,label)

  !> Checks that two arrays are of the same size.

    use glimmer_log

    real(dp),dimension(:,:),intent(in) :: array1 !< The array 1 to be checked
    real(dp),dimension(:,:),intent(in) :: array2 !< The array 2 to be checked
    character(*),intent(in),optional :: label    !< Optional label, to facilitate bug tracking if the check fails.

    if ((size(array1,1)/=size(array2,1)).or.(size(array1,2)/=size(array2,2))) then
      if (present(label)) then
        call write_log('Non-conformal arrays. Label: '//label,GM_FATAL,__FILE__,__LINE__)
      else
        call write_log('ERROR: Non-conformal arrays. No label',GM_FATAL,__FILE__,__LINE__)
      endif
    endif

  end subroutine check_conformal_2d_real

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> compute horizontal sum for each vertical level
  !!
  !! Calculates the sum of a given three-dimensional field at each
  !! level. The vertical coordinate of the input is the first index of
  !! the array.
  !! \return
  !! A one-dimensional array of the same size as the first dimension of
  !! inp is returned, containing the sum of inp for 
  !! each level.
  function hsum(inp)


    implicit none

    real(dp),dimension(:,:,:),intent(in) :: inp !< The input array. The first index is the vertical, the othe two horizontal.
    real(dp),dimension(size(inp,dim=1))  :: hsum
  
    integer up

    do up=1,size(inp,dim=1)
       hsum(up) = sum(inp(up,:,:))
    end do

  end function hsum

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Calculates the sum of a given two-dimensional field along one axis.
  !! Within GLIMMER, this function calculates the mean vertical profile
  !! in a 2D vertical slice. 
  !! \return
  !! A one-dimensional array of the same size as the first dimension of
  !! inp is returned, containing the sum of inp for 
  !! each row.

  function lsum(inp)


    implicit none

    real(dp),dimension(:,:), intent(in) :: inp !< Input array
    real(dp),dimension(size(inp,dim=1)) :: lsum
    
    lsum = sum(inp(:,:),dim=2)

  end function lsum

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Tridiagonal solver. All input/output arrays should have the 
  !! same number of elements.

  subroutine tridiag(a,b,c,x,y)

    real(dp), dimension(:), intent(in)  :: a !< Lower diagonal; a(1) is ignored.
    real(dp), dimension(:), intent(in)  :: b !< Center diagonal
    real(dp), dimension(:), intent(in)  :: c !< Upper diagonal; c(n) is ignored.
    real(dp), dimension(:), intent(out) :: x !< Unknown vector
    real(dp), dimension(:), intent(in)  :: y !< Right-hand side

    real(dp), dimension(size(a)) :: aa
    real(dp), dimension(size(a)) :: bb

    integer :: n,i

    n = size(a)

    aa(1) = c(1)/b(1)
    bb(1) = y(1)/b(1)

    do i = 2,n
       aa(i) = c(i) / (b(i)-a(i)*aa(i-1))
       bb(i) = (y(i)-a(i)*bb(i-1)) / (b(i)-a(i)*aa(i-1))
    end do

    x(n) = bb(n)

    do i = n-1,1,-1
       x(i) = bb(i)-aa(i)*x(i+1)
    end do

  end subroutine tridiag

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Strip leading and trailing single or double quotes from a string.
  !! Returns the string with these quote characters removed, or a copy of the original
  !! string if it didn't have leading and trailing quotes.
  function strip_quotes(str) result(stripped_str)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: stripped_str

    integer :: n

    n = len_trim(str)
    if ((str(1:1) == "'" .and. str(n:n) == "'") .or. &
         (str(1:1) == '"' .and. str(n:n) == '"')) then
       ! This is a quoted string; return the string with the leading and trailing quotes removed
       stripped_str = str(2:(n-1))
    else
       stripped_str = str
    end if
  end function strip_quotes

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_lsrf_usrf(thck, topg, eus, lsrf, usrf)

    ! Calculate the elevation of the lower and upper surface of the ice,
    ! given the thickness and bed topography.

    ! Note: This subroutine computes over all grid cells, not just locally owned.
    !       Output will be correct in halos only if the input is correct.
    !       Generally the units will be meters, but the output will be correct
    !        as long as the units are mutually consistent.

    use glimmer_physcon, only : rhoi, rhoo

    implicit none

    real(dp), intent(in),  dimension(:,:) :: thck   !> ice thickness
    real(dp), intent(in),  dimension(:,:) :: topg   !> bedrock topography elevation
    real(dp), intent(in)                  :: eus    !> global sea level

    real(dp), intent(out), dimension(:,:) :: lsrf   !> lower ice surface elevation
    real(dp), intent(out), dimension(:,:) :: usrf   !> upper ice surface elevation

    ! Compute lsrf by considering whether the ice is floating or not
    ! For ice-free land, lsrf = topg
    ! For ice-free ocean, lsrf = 0

    where (topg - eus < (-rhoi/rhoo) * thck)
       lsrf = eus - (rhoi/rhoo) * thck
    elsewhere
       lsrf = topg
    end where

    ! Compute usrf
    usrf = lsrf + thck

  end subroutine calc_lsrf_usrf

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine point_diag_integer_2d(&
       field,  field_string,   &
       ipt,    jpt,    rpt,    &
       irange, jrange,         &
       fmt_str)

    use cism_parallel, only: this_rank

    ! Write the value of a field in a small neighborhood around a central point

    ! input/output arguments

    integer, dimension(:,:), intent(in) :: &
         field                        ! field to be written to stdout

    character(*), intent(in) :: &
         field_string                 ! field name or a similar string, e.g. 'ice_mask'

    integer, intent(in) :: &
         ipt, jpt, rpt                ! i and j coordinates and processor rank of central point

    integer, intent(in) :: &
         irange,                    & ! range of diagnostic output in i direction
         jrange                       ! range of diagnostic output in j direction

    character(len=*), intent(in), optional :: &
         fmt_str                      ! string for fortran formatted output, default = '(f10.3)'

    ! local variables

    integer :: i, j
    integer :: ilo, jlo, ihi, jhi
    integer :: left_range, right_range
    character(8) :: fmtstr               ! string for fortran formatted output

    if (this_rank == rpt) then

       if (present(fmt_str)) then
          fmtstr = trim(fmt_str)
!          write(iulog,*) 'In point_diag, fmtstr = ', trim(fmtstr)
       else
          fmtstr = '(i10)'
!          write(iulog,*) 'In point_diag, fmtstr = ', trim(fmtstr)
       endif

       if (mod(irange,2) == 0) then  ! even number
          ilo = ipt - irange/2
          ihi = ipt + irange/2 - 1
       else   ! odd number
          ilo = ipt - (irange-1)/2
          ihi = ipt + (irange-1)/2
       endif

       if (mod(jrange,2) == 0) then  ! even number
          jlo = jpt - jrange/2
          jhi = jpt + jrange/2 - 1
       else   ! odd number
          jlo = jpt - (jrange-1)/2
          jhi = jpt + (jrange-1)/2
       endif

       ! Make sure all points are in bounds
       ilo = max(ilo,1)
       ihi = min(ihi,size(field,1))
       jlo = max(jlo,1)
       jhi = min(jhi,size(field,2))

       ! Write to stdout

       write(iulog,*) ' '
       write(iulog,*) trim(field_string), ': i, j, rank =', ipt, jpt, rpt
       do j = jhi, jlo, -1    ! top to bottom
          do i = ilo, ihi     ! left to right
             write(iulog,fmtstr,advance='no') field(i,j)
          enddo
          write(iulog,*) ' '
       enddo

    endif  ! this_rank = rtest

  end subroutine point_diag_integer_2d

!--------------------------------------------------------------------------

  subroutine point_diag_logical_2d(&
       field,  field_string,   &
       ipt,    jpt,    rpt,    &
       irange, jrange,         &
       fmt_str)

    use cism_parallel, only: this_rank

    ! Write the value of a field in a small neighborhood around a central point

    ! input/output arguments

    logical, dimension(:,:), intent(in) :: &
         field                        ! field to be written to stdout

    character(*), intent(in) :: &
         field_string                 ! field name or a similar string, e.g. 'ice_mask'

    integer, intent(in) :: &
         ipt, jpt, rpt                ! i and j coordinates and processor rank of central point

    integer, intent(in) :: &
         irange,                    & ! range of diagnostic output in i direction
         jrange                       ! range of diagnostic output in j direction

    character(len=*), intent(in), optional :: &
         fmt_str                      ! string for fortran formatted output, default = '(f10.3)'

    ! local variables

    integer :: i, j
    integer :: ilo, jlo, ihi, jhi
    integer :: left_range, right_range
    character(8) :: fmtstr               ! string for fortran formatted output

    if (this_rank == rpt) then

       if (present(fmt_str)) then
          fmtstr = trim(fmt_str)
!          write(iulog,*) 'In point_diag, fmtstr = ', trim(fmtstr)
       else
          fmtstr = '(L10)'
!          write(iulog,*) 'In point_diag, fmtstr = ', trim(fmtstr)
       endif

       if (mod(irange,2) == 0) then  ! even number
          ilo = ipt - irange/2
          ihi = ipt + irange/2 - 1
       else   ! odd number
          ilo = ipt - (irange-1)/2
          ihi = ipt + (irange-1)/2
       endif

       if (mod(jrange,2) == 0) then  ! even number
          jlo = jpt - jrange/2
          jhi = jpt + jrange/2 - 1
       else   ! odd number
          jlo = jpt - (jrange-1)/2
          jhi = jpt + (jrange-1)/2
       endif

       ! Make sure all points are in bounds
       ilo = max(ilo,1)
       ihi = min(ihi,size(field,1))
       jlo = max(jlo,1)
       jhi = min(jhi,size(field,2))

       ! Write to stdout

       write(iulog,*) ' '
       write(iulog,*) trim(field_string), ': i, j, rank =', ipt, jpt, rpt
       do j = jhi, jlo, -1    ! top to bottom
          do i = ilo, ihi     ! left to right
             write(iulog,fmtstr,advance='no') field(i,j)
          enddo
          write(iulog,*) ' '
       enddo

    endif  ! this_rank = rtest

  end subroutine point_diag_logical_2d

!--------------------------------------------------------------------------

  subroutine point_diag_real8_2d(&
       field,  field_string,   &
       ipt,    jpt,    rpt,    &
       irange, jrange,         &
       fmt_str)

    use cism_parallel, only: this_rank

    ! Write the value of a field in a small neighborhood around a central point

    ! input/output arguments

    real(dp), dimension(:,:), intent(in) :: &
         field                        ! field to be written to stdout

    character(*), intent(in) :: &
         field_string                 ! field name or a similar string, e.g. 'thck' or 'thck (m)'

    integer, intent(in) :: &
         ipt, jpt, rpt                ! i and j coordinates and processor rank of central point

    integer, intent(in) :: &
         irange,                    & ! range of diagnostic output in i direction
         jrange                       ! range of diagnostic output in j direction

    character(len=*), intent(in), optional :: &
         fmt_str                      ! string for fortran formatted output, default = '(f10.3)'

    ! local variables

    integer :: i, j
    integer :: ilo, jlo, ihi, jhi
    character(8) :: fmtstr               ! string for fortran formatted output

    if (this_rank == rpt) then

       if (present(fmt_str)) then
          fmtstr = trim(fmt_str)
       else
          fmtstr = '(f10.3)'
       endif

!       write(iulog,*) 'In point_diag, fmtstr = ', trim(fmtstr)

       if (mod(irange,2) == 0) then  ! even number
          ilo = ipt - irange/2
          ihi = ipt + irange/2 - 1
       else   ! odd number
          ilo = ipt - (irange-1)/2
          ihi = ipt + (irange-1)/2
       endif

       if (mod(jrange,2) == 0) then  ! even number
          jlo = jpt - jrange/2
          jhi = jpt + jrange/2 - 1
       else   ! odd number
          jlo = jpt - (jrange-1)/2
          jhi = jpt + (jrange-1)/2
       endif

       ! Make sure all points are in bounds
       ilo = max(ilo,1)
       ihi = min(ihi,size(field,1))
       jlo = max(jlo,1)
       jhi = min(jhi,size(field,2))

       ! Write to stdout
       write(iulog,*) ' '
       write(iulog,*) trim(field_string), ': i, j, rank =', ipt, jpt, rpt
       do j = jhi, jlo, -1    ! top to bottom
          do i = ilo, ihi     ! left to right
             write(iulog,fmtstr,advance='no') field(i,j)
          enddo
          write(iulog,*) ' '
       enddo

    endif  ! this_rank = rpt

  end subroutine point_diag_real8_2d

!--------------------------------------------------------------------------

  subroutine double_to_binary(&
       x, binary_str, binary_full, binary_sign, binary_exponent, binary_mantissa)

    ! Find the internal binary representation of a double-precision floating point number
    ! Based on the IEEE-754 standard

    use glimmer_global, only: dp, i8
    implicit none

    real(dp), intent(in) :: x
    character(len=64), intent(out) :: binary_str    ! string representation of the binary number

    integer(i8), intent(out), optional :: binary_full      ! 64 bits
    integer, intent(out), optional :: binary_sign          ! 1 bit
    integer, intent(out), optional :: binary_exponent      ! 11 bits
    integer(i8), intent(out), optional :: binary_mantissa  ! 52 bits

    integer :: i
    character(len=1) :: bin(64)
    integer (i8) :: binary_number
    integer :: sign_bit
    integer :: exponent_bits
    integer :: mantissa_bits

    logical :: verbose_binary = .false.

    ! Transfer the double value into a 64-bit integer
    binary_number = transfer(x, binary_number)

    ! Get the sign bit (bit 64)
    sign_bit = ishft(binary_number, -63) .and. 1

    ! Get the exponent bits (bits 63–53)
    exponent_bits = ishft(binary_number, -52) .and. Z'7FF'

    ! Extract mantissa (fraction) bits (bits 52–1)
    mantissa_bits = binary_number .and. Z'FFFFFFFFFFFFF'

    if (present(binary_full)) binary_full = binary_number
    if (present(binary_sign)) binary_sign = sign_bit
    if (present(binary_exponent)) binary_exponent = exponent_bits
    if (present(binary_mantissa)) binary_mantissa = mantissa_bits

    if (verbose_binary) then
       write(iulog,*) ' '
       write(iulog,*) 'x =', x
       write(iulog,*) 'IEEE-754 double precision representation of x:'
       write(iulog,*) 'Sign bit: ', sign_bit
       write(iulog,*) 'Exponent (11 bits):', exponent_bits
       write(iulog,*) 'Mantissa (52 bits):', mantissa_bits
    endif

    ! Convert full 64-bit integer to a binary string
    do i = 1, 64
        if (btest(binary_number, 64 - i)) then
            bin(i) = '1'
        else
            bin(i) = '0'
        end if
    end do

    binary_str = concat(bin)
    if (verbose_binary) then
       write(iulog,*) 'Full 64-bit binary:'
       write(iulog,*), ' ', binary_str
    endif

  end subroutine double_to_binary

!--------------------------------------------------------------------------

  subroutine cubic_solver(&
       a, b, c, d, &
       x1,         &
       x2_r, x2_i, &
       x3_r, x3_i)

    !------------------------------------------------
    !Note: This subroutine is not currently called, but is here in case
    !      it turns out to be useful.
    !TODO - Pass 3 complex roots in and out?
    !------------------------------------------------
    ! Find the real root of a cubic equation:
    !
    !    ax^3 + bx^2 + cx = d = 0
    !
    ! Do this by making the substitution
    !
    !    x = y - b/(3a)
    !
    ! to convert to a depressed cubic:
    !
    !    y^3 + py + q = 0
    !
    ! where p = (1/a) * (c - b^2/(3a))
    !       q = (1/a) * (d + 2b^3/(27a^2) - bc/(3a))
    !
    !------------------------------------------------

    use glimmer_physcon, only: pi

    real(dp), intent(in) ::  &
         a, b, c, d       ! coefficients of cubic equation
                          ! assumed to be real

    real(dp), intent(out) ::  &
         x1               ! real solution of cubic equation

    real(dp), intent(out), optional ::  &
         x2_r, x2_i,    & ! other solutions of cubic equation
         x3_r, x3_i       ! could be either real or complex

    real(dp) :: &
         p, q             ! coefficients of depressed cubic

    real(dp) :: &
         Delta            ! discriminant

    real(dp) :: &
         y1,            & ! solutions of depressed cubic
         y2_r, y2_i,    & !
         y3_r, y3_i

    real(dp) :: &
         u, v,          & ! some intermediate factors
         fu, fv,        &
         phi

    real(dp), parameter :: &
         p333 = 1.d0/3.d0

    logical, parameter :: verbose_cubic = .false.

    ! compute coefficients of depressed cubic, y^3 + py + q = 0

    p = (3.d0*c/a - (b/a)**2) / 3.d0
    q = (2.d0*(b/a)**3 - 9.d0*b*c/(a*a) + 27.d0*d/a) / 27.d0

    ! compute the discriminant
    Delta = (p/3.d0)**3 + (q/2.d0)**2

    if (verbose_cubic) then
       write(iulog,*) 'Delta =', Delta
       if (Delta > 0.d0) then
          write(iulog,*) 'One real root, 2 complex conjugate'
       elseif (Delta == 0.d0) then
          write(iulog,*) 'Three real roots of which at least two are equal'
       elseif (Delta < 0.d0) then
          write(iulog,*) 'Three distinct real roots'
       endif
    endif

    if (Delta >= 0.d0) then

       if (Delta > 0.d0) then    ! one real root, two complex roots
          fu = -q/2.d0 + sqrt(Delta)
          fv = -q/2.d0 - sqrt(Delta)
       else  ! Delta = 0; three real roots of which at least two are equal
          fu = -q/2.d0
          fv = fu
       endif

       ! some logic to avoid taking cube roots of negative numbers
       if (fu >= 0.d0) then
          u = fu**p333
       else
          u = -(-fu)**p333
       endif

       if (fv >= 0.d0) then
          v = fv**p333
       else
          v = -(-fv)**p333
       endif

       ! form solutions of depressed cubic
       y1 = u + v       ! real
       y2_r = -(u+v)/2.d0
       y2_i =  (u-v)*sqrt(3.d0)/2.d0
       y3_r = -(u+v)/2.d0
       y3_r = -(u-v)*sqrt(3.d0)/2.d0

       if (verbose_cubic) then
          write(iulog,*) 'a, b, c, d:', a, b, c, d
          write(iulog,*) 'p, q:', p, q
          write(iulog,*) 'y1 =', y1
          write(iulog,*) 'x1 =', x1
       endif

    else  ! Delta < 0; three distinct real roots
          ! use a trigonometric formulation

       phi = acos(-q/(2.d0*sqrt(abs(p)**3/27.d0)))

       y1 =    2.d0 * sqrt(abs(p)/3.d0) * cos(phi/3.d0)
       y2_r = -2.d0 * sqrt(abs(p)/3.d0) * cos((phi+pi)/3.d0)
       y2_i =  0.d0
       y3_r = -2.d0 * sqrt(abs(p)/3.d0) * cos((phi-pi)/3.d0)
       y3_i =  0.d0

       if (verbose_cubic) then
          write(iulog,*) 'a, b, c, d:', a, b, c, d
          write(iulog,*) 'p, q:', p, q
          write(iulog,*) 'y1, y2, y3 =', y1, y2_r, y3_r
          write(iulog,*) 'b/3a =', b/(3.d0*a)
          write(iulog,*) 'x1 =', y1 - b/(3.d0*a)
          write(iulog,*) 'x2 =', y2_r - b/(3.d0*a)
          write(iulog,*) 'x3 =', y3_r - b/(3.d0*a)
       endif

    endif

    ! Recover the solutions
    ! Mostly likely we are only interested in x1, but compute the others if requested

    x1 = y1 - b/(3.d0*a)

    if (present(x2_r) .and. present(x2_i) .and. present(x3_r) .and. present(x3_i)) then
       x2_r = y2_r - b/(3.d0*a)
       x2_i = y2_i
       x3_r = y3_r - b/(3.d0*a)
       x3_i = y3_i
    endif

  end subroutine cubic_solver

!--------------------------------------------------------------------------

  pure function concat(arr) result(str)
    ! Turn a character array into a string

    character(len=*), intent(in) :: arr(:)
    character(len=size(arr)) :: str
    integer :: k

    do k = 1, size(arr)
       str(k:k) = arr(k)
    end do

  end function concat

!****************************************************************************

end module glimmer_utils

!**************************************************************************** 
