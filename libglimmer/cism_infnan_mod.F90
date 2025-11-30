! Flag representing compiler support of Fortran 2003's
! ieee_arithmetic intrinsic module.
#if defined CPRIBM || defined CPRPGI || defined CPRINTEL ||  defined CPRCRAY || defined CPRNAG
#define HAVE_IEEE_ARITHMETIC
#endif

!WHL, Nov. 2025: Adapted from shr_infnan_mod.F90, part of CESM shared code
!                Changed 'shr' to 'cism' to avoid name conflicts with shared code
!                I kept only the r8 interfaces (not r4, i8 or i4).
!                Also, I assumed that all input arrays are 1d.

module cism_infnan_mod
!---------------------------------------------------------------------
! Module to test for IEEE Inf and NaN values, which also provides a
! method of setting +/-Inf and signaling or quiet NaN.
!
! All functions are elemental, and thus work on arrays.
!---------------------------------------------------------------------
! To test for these values, just call the corresponding function, e.g:
!
!   var_is_nan = cism_infnan_isnan(x)
!
! You can also use it on arrays:
!
!   array_contains_nan = any(cism_infnan_isnan(my_array))
!
!---------------------------------------------------------------------
! To generate these values, assign one of the provided derived-type
! variables to a real:
!
!   use cism_infnan_mod, only: nan => cism_infnan_nan, &
!                             inf => cism_infnan_inf, &
!                             assignment(=)
!   real(r4) :: my_nan
!   real(r8) :: my_inf_array(2,2)
!   my_nan = nan
!   my_inf_array = inf
!
! Keep in mind that "cism_infnan_nan" and "cism_infnan_inf" cannot be
! passed to functions that expect real arguments. To pass a real
! NaN, you will have to use cism_infnan_nan to set a local real of
! the correct kind.
!---------------------------------------------------------------------

  use glimmer_global, only: r4 => sp, r8 => dp
  use glimmer_global, only: i4, i8
!!use shr_kind_mod, only: &
!!     r4 => SHR_KIND_R4, &
!!     r8 => SHR_KIND_R8

#ifdef HAVE_IEEE_ARITHMETIC

! If we have IEEE_ARITHMETIC, the NaN test is provided for us.
use, intrinsic :: ieee_arithmetic, only: &
     cism_infnan_isnan => ieee_is_nan

#else

! Integers of correct size for bit patterns below.
!!use shr_kind_mod, only: i4 => shr_kind_i4, i8 => shr_kind_i8

#endif

implicit none
private
save

! Test functions for NaN/Inf values.
public :: cism_infnan_isnan
public :: cism_infnan_isinf
public :: cism_infnan_isposinf
public :: cism_infnan_isneginf

! Locally defined isnan.
#ifndef HAVE_IEEE_ARITHMETIC
interface cism_infnan_isnan
   ! TYPE double,real
   module procedure cism_infnan_isnan_r8
end interface
#endif

interface cism_infnan_isinf
   ! TYPE double,real
   module procedure cism_infnan_isinf_r8
end interface

interface cism_infnan_isposinf
   ! TYPE double,real
   module procedure cism_infnan_isposinf_r8
end interface

interface cism_infnan_isneginf
   ! TYPE double,real
   module procedure cism_infnan_isneginf_r8
end interface

! Derived types for generation of NaN/Inf
! Even though there's no reason to "use" the types directly, some compilers
! might have trouble with an object being used without its type.
public :: cism_infnan_nan_type
public :: cism_infnan_inf_type
public :: assignment(=)
public :: cism_infnan_to_r4
public :: cism_infnan_to_r8

! Type representing Not A Number.
type :: cism_infnan_nan_type
   logical :: quiet = .false.
end type cism_infnan_nan_type

! Type representing +/-Infinity.
type :: cism_infnan_inf_type
   logical :: positive = .true.
end type cism_infnan_inf_type

! Allow assigning reals to NaN or Inf.
interface assignment(=)
   ! TYPE double,real
   ! DIMS 0,1,2,3,4,5,6,7
   module procedure set_nan_1d_r8
   ! TYPE double,real
   ! DIMS 0,1,2,3,4,5,6,7
   module procedure set_inf_1d_r8
end interface

! Conversion functions.
interface cism_infnan_to_r8
   module procedure nan_r8
   module procedure inf_r8
end interface

interface cism_infnan_to_r4
   module procedure nan_r4
   module procedure inf_r4
end interface

! Initialize objects of NaN/Inf type for other modules to use.

! Default NaN is signaling, but also provide snan and qnan to choose
! explicitly.
type(cism_infnan_nan_type), public, parameter :: cism_infnan_nan = &
     cism_infnan_nan_type(.false.)
type(cism_infnan_nan_type), public, parameter :: cism_infnan_snan = &
     cism_infnan_nan_type(.false.)
type(cism_infnan_nan_type), public, parameter :: cism_infnan_qnan = &
     cism_infnan_nan_type(.true.)

! Default Inf is positive, but provide posinf to go with neginf.
type(cism_infnan_inf_type), public, parameter :: cism_infnan_inf = &
     cism_infnan_inf_type(.true.)
type(cism_infnan_inf_type), public, parameter :: cism_infnan_posinf = &
     cism_infnan_inf_type(.true.)
type(cism_infnan_inf_type), public, parameter :: cism_infnan_neginf = &
     cism_infnan_inf_type(.false.)

! Bit patterns for implementation without ieee_arithmetic.
! Note that in order to satisfy gfortran's range check, we have to use
! ibset to set the sign bit from a BOZ pattern.
#ifndef HAVE_IEEE_ARITHMETIC
! Single precision.
integer(i4), parameter :: ssnan_pat = int(Z'7FA00000',i4)
integer(i4), parameter :: sqnan_pat = int(Z'7FC00000',i4)
integer(i4), parameter :: sposinf_pat = int(Z'7F800000',i4)
integer(i4), parameter :: sneginf_pat = ibset(sposinf_pat,bit_size(1_i4)-1)
! Double precision.
integer(i8), parameter :: dsnan_pat = int(Z'7FF4000000000000',i8)
integer(i8), parameter :: dqnan_pat = int(Z'7FF8000000000000',i8)
integer(i8), parameter :: dposinf_pat = int(Z'7FF0000000000000',i8)
integer(i8), parameter :: dneginf_pat = ibset(dposinf_pat,bit_size(1_i8)-1)
#endif

contains

!---------------------------------------------------------------------
! TEST FUNCTIONS
!---------------------------------------------------------------------
! The "isinf" function simply calls "isposinf" and "isneginf".
!---------------------------------------------------------------------

! TYPE double,real
elemental function cism_infnan_isinf_r8(x) result(isinf)
  real(r8), intent(in) :: x
  logical :: isinf

  isinf = cism_infnan_isposinf(x) .or. cism_infnan_isneginf(x)

end function cism_infnan_isinf_r8

#ifdef HAVE_IEEE_ARITHMETIC

!---------------------------------------------------------------------
! The "isposinf" and "isneginf" functions get the IEEE class of a
! real, and test to see if the class is equal to ieee_positive_inf
! or ieee_negative_inf.
!---------------------------------------------------------------------

! TYPE double,real
elemental function cism_infnan_isposinf_r8(x) result(isposinf)
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_class, &
       ieee_positive_inf, &
       operator(==)
  real(r8), intent(in) :: x
  logical :: isposinf

  isposinf = (ieee_positive_inf == ieee_class(x))

end function cism_infnan_isposinf_r8

! TYPE double,real
elemental function cism_infnan_isneginf_r8(x) result(isneginf)
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_class, &
       ieee_negative_inf, &
       operator(==)
  real(r8), intent(in) :: x
  logical :: isneginf

  isneginf = (ieee_negative_inf == ieee_class(x))

end function cism_infnan_isneginf_r8

#else
! Don't have ieee_arithmetic.

!!#ifdef CPRGNU   !WHL - Assume this is true
! NaN testing on gfortran.
! TYPE double,real
elemental function cism_infnan_isnan_r8(x) result(is_nan)
  real(r8), intent(in) :: x
  logical :: is_nan

  is_nan = isnan(x)

end function cism_infnan_isnan_r8
! End GNU section.
!!#endif

!---------------------------------------------------------------------
! The "isposinf" and "isneginf" functions just test against a known
! bit pattern if we don't have ieee_arithmetic.
!---------------------------------------------------------------------

! TYPE double,real
elemental function cism_infnan_isposinf_r8(x) result(isposinf)
  real(r8), intent(in) :: x
  logical :: isposinf
!!#if ({ITYPE} == TYPEREAL)
!!  integer(i4), parameter :: posinf_pat = sposinf_pat
!!#else
  integer(i8), parameter :: posinf_pat = dposinf_pat
!!#endif

  isposinf = (x == transfer(posinf_pat,x))

end function cism_infnan_isposinf_r8

! TYPE double,real
elemental function cism_infnan_isneginf_r8(x) result(isneginf)
  real(r8), intent(in) :: x
  logical :: isneginf
!!#if ({ITYPE} == TYPEREAL)
!!  integer(i4), parameter :: neginf_pat = sneginf_pat
!!#else
  integer(i8), parameter :: neginf_pat = dneginf_pat
!!#endif

  isneginf = (x == transfer(neginf_pat,x))

end function cism_infnan_isneginf_r8

! End ieee_arithmetic conditional.
#endif

!---------------------------------------------------------------------
! GENERATION FUNCTIONS
!---------------------------------------------------------------------
! Two approaches for generation of NaN and Inf values:
!   1. With Fortran 2003, use the ieee_value intrinsic to get a value
!      from the corresponding class. These are:
!       - ieee_signaling_nan
!       - ieee_quiet_nan
!       - ieee_positive_inf
!       - ieee_negative_inf
!   2. Without Fortran 2003, set the IEEE bit patterns directly.
!      Use BOZ literals to get an integer with the correct bit
!      pattern, then use "transfer" to transfer those bits into a
!      real.
!---------------------------------------------------------------------

! TYPE double,real
! DIMS 0,1,2,3,4,5,6,7
pure subroutine set_nan_1d_r8(output, nan)
#ifdef HAVE_IEEE_ARITHMETIC
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_signaling_nan, &
       ieee_quiet_nan, &
       ieee_value
#else
!!#if ({ITYPE} == TYPEREAL)
!!  integer(i4), parameter :: snan_pat = ssnan_pat
!!  integer(i4), parameter :: qnan_pat = sqnan_pat
!!#else
  integer(i8), parameter :: snan_pat = dsnan_pat
  integer(i8), parameter :: qnan_pat = dqnan_pat
!!#endif
#endif
!!  real(r8), intent(out) :: output{DIMSTR}
  real(r8), intent(out) :: output
  type(cism_infnan_nan_type), intent(in) :: nan

  ! Use scalar temporary for performance reasons, to reduce the cost of
  ! the ieee_value call.
  real(r8) :: tmp

#ifdef HAVE_IEEE_ARITHMETIC
  if (nan%quiet) then
     tmp = ieee_value(tmp, ieee_quiet_nan)
  else
     tmp = ieee_value(tmp, ieee_signaling_nan)
  end if
#else
  if (nan%quiet) then
     tmp = transfer(qnan_pat, tmp)
  else
     tmp = transfer(snan_pat, tmp)
  end if
#endif

  output = tmp

end subroutine set_nan_1d_r8

! TYPE double,real
! DIMS 0,1,2,3,4,5,6,7
pure subroutine set_inf_1d_r8(output, inf)
#ifdef HAVE_IEEE_ARITHMETIC
  use, intrinsic :: ieee_arithmetic, only: &
       ieee_positive_inf, &
       ieee_negative_inf, &
       ieee_value
#else
!!#if ({ITYPE} == TYPEREAL)
!!  integer(i4), parameter :: posinf_pat = sposinf_pat
!!  integer(i4), parameter :: neginf_pat = sneginf_pat
!!#else
  integer(i8), parameter :: posinf_pat = dposinf_pat
  integer(i8), parameter :: neginf_pat = dneginf_pat
!!#endif
#endif
!!  real(r8), intent(out) :: output{DIMSTR}
  real(r8), intent(out) :: output
  type(cism_infnan_inf_type), intent(in) :: inf

  ! Use scalar temporary for performance reasons, to reduce the cost of
  ! the ieee_value call.
  real(r8) :: tmp

#ifdef HAVE_IEEE_ARITHMETIC
  if (inf%positive) then
     tmp = ieee_value(tmp,ieee_positive_inf)
  else
     tmp = ieee_value(tmp,ieee_negative_inf)
  end if
#else
  if (inf%positive) then
     tmp = transfer(posinf_pat, tmp)
  else
     tmp = transfer(neginf_pat, tmp)
  end if
#endif

  output = tmp

end subroutine set_inf_1d_r8

!---------------------------------------------------------------------
! CONVERSION INTERFACES.
!---------------------------------------------------------------------
! Function methods to get reals from nan/inf types.
!---------------------------------------------------------------------

pure function nan_r8(nan) result(output)
  class(cism_infnan_nan_type), intent(in) :: nan
  real(r8) :: output

!!  output = nan
  !WHL kluge
  output = 0._r8

end function nan_r8

pure function nan_r4(nan) result(output)
  class(cism_infnan_nan_type), intent(in) :: nan
  real(r4) :: output

!!  output = nan
  !WHL kluge
  output = 0._r8

end function nan_r4

pure function inf_r8(inf) result(output)
  class(cism_infnan_inf_type), intent(in) :: inf
  real(r8) :: output

!!  output = inf
  !WHL kluge
  output = 0._r8

end function inf_r8

pure function inf_r4(inf) result(output)
  class(cism_infnan_inf_type), intent(in) :: inf
  real(r4) :: output

!!  output = inf
  !WHL kluge
  output = 0._r8

end function inf_r4

end module cism_infnan_mod
