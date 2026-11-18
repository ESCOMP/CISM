!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   glissade_log_config.F90 - part of the Community Ice Sheet Model (CISM)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glissade_log_config

  ! Helpers that write the model configuration to the log file in the style of a config file:
  !
  !     [velocity]
  !       stokes_approx                  = 4             # Depth-integrated viscosity (glissade_velo_higher)
  !       max_slope                      = 1.0
  !       effstrain_min                  = 1.0e-08       # yr^-1
  !
  ! The print_<section> subroutines in glissade_setup.F90 call these helpers, so that every section
  ! is written in the same format. Conventions:
  !   - log_section writes a blank line and the section header.
  !   - One line per key: the key name exactly as in the config file, ' = ', the value, and an
  !     optional comment after '#' (the description of an option value, the units of a parameter,
  !     or a note, e.g. when the code has changed a setting).
  !   - Options given by name are written by name, with the description of the value as the
  !     comment (log_option with a choices array; see option_choice in glimmer_config).
  !   - Integers are written without padding. Logicals are written as .true. or .false.
  !   - Reals are written by real_to_string (see below): at most 12 significant digits, without
  !     trailing zeros; fixed notation for 1.e-6 <= |x| < 1.e6, scientific notation otherwise.
  !   - Information that is not a config key (e.g. a derived quantity or an explanation) is written
  !     as a comment line with log_comment.

  use glimmer_global, only: dp
  use glimmer_log, only: write_log
  use glimmer_config, only: option_choice, option_name, option_description

  implicit none
  private

  public :: log_section, log_comment, log_option, log_param, log_flag, log_string, log_list, real_to_string

  interface log_param
     module procedure log_param_real, log_param_int
  end interface log_param

  interface log_option
     module procedure log_option_choice, log_option_int
  end interface log_option

  integer, parameter :: indent_key = 4       ! spaces before a key
  integer, parameter :: key_width = 30       ! keys are padded to this width (longer keys are not)
  integer, parameter :: value_width = 13     ! values followed by a comment are padded to this width
  integer, parameter :: default_digits = 12  ! default number of significant digits for reals
  real(dp), parameter :: fixed_min = 1.0d-6  ! smallest magnitude written in fixed notation
  real(dp), parameter :: fixed_max = 1.0d6   ! magnitudes >= fixed_max are written in scientific notation

contains

!--------------------------------------------------------------------------------

  subroutine log_section(name)

    ! Write a blank line and a section header, e.g. '  [velocity]'

    character(len=*), intent(in) :: name

    call write_log(' ')
    call write_log('  [' // trim(name) // ']')

  end subroutine log_section

!--------------------------------------------------------------------------------

  subroutine log_comment(text)

    ! Write a comment line inside a section, e.g. '    # Icebergs will be removed'

    character(len=*), intent(in) :: text

    call write_log(repeat(' ', indent_key) // '# ' // trim(text))

  end subroutine log_comment

!--------------------------------------------------------------------------------

  subroutine log_option_choice(key, value, choices, note)

    ! Write an option given by name: the name of its value, with the description as the comment.
    ! An optional note is added to the comment. A value that is not one of the choices (which
    ! the config reader does not allow) is written as an integer.

    character(len=*), intent(in) :: key
    integer, intent(in) :: value
    type(option_choice), dimension(:), intent(in) :: choices
    character(len=*), intent(in), optional :: note

    character(len=32) :: name
    character(len=256) :: comment

    name = option_name(choices, value)
    if (len_trim(name) == 0) then
       call write_key_value(key, int_to_string(value), 'not a valid value')
       return
    endif
    comment = option_description(choices, value)
    if (present(note)) then
       if (len_trim(comment) > 0) then
          comment = trim(comment) // '; ' // note
       else
          comment = note
       endif
    endif
    if (len_trim(comment) > 0) then
       call write_key_value(key, trim(name), trim(comment))
    else
       call write_key_value(key, trim(name))
    endif

  end subroutine log_option_choice

!--------------------------------------------------------------------------------

  subroutine log_option_int(key, value, description)

    ! Write an integer option and, optionally, the description of its value

    character(len=*), intent(in) :: key
    integer, intent(in) :: value
    character(len=*), intent(in), optional :: description

    call write_key_value(key, int_to_string(value), description)

  end subroutine log_option_int

!--------------------------------------------------------------------------------

  subroutine log_param_real(key, value, comment, digits, style)

    ! Write a real parameter; see real_to_string for digits and style

    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value
    character(len=*), intent(in), optional :: comment    ! e.g. the units
    integer, intent(in), optional :: digits
    character(len=*), intent(in), optional :: style

    call write_key_value(key, real_to_string(value, digits, style), comment)

  end subroutine log_param_real

!--------------------------------------------------------------------------------

  subroutine log_param_int(key, value, comment)

    ! Write an integer parameter (a count or an index rather than an option)

    character(len=*), intent(in) :: key
    integer, intent(in) :: value
    character(len=*), intent(in), optional :: comment

    call write_key_value(key, int_to_string(value), comment)

  end subroutine log_param_int

!--------------------------------------------------------------------------------

  subroutine log_flag(key, value, comment)

    ! Write a logical option as .true. or .false.

    character(len=*), intent(in) :: key
    logical, intent(in) :: value
    character(len=*), intent(in), optional :: comment

    if (value) then
       call write_key_value(key, '.true.', comment)
    else
       call write_key_value(key, '.false.', comment)
    endif

  end subroutine log_flag

!--------------------------------------------------------------------------------

  subroutine log_string(key, value, comment)

    ! Write a character value, e.g. a file name

    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value
    character(len=*), intent(in), optional :: comment

    call write_key_value(key, trim(value), comment)

  end subroutine log_string

!--------------------------------------------------------------------------------

  subroutine log_list(key, values, comment, digits)

    ! Write a list of reals separated by spaces, as in the config file (e.g. sigma_levels)

    character(len=*), intent(in) :: key
    real(dp), dimension(:), intent(in) :: values
    character(len=*), intent(in), optional :: comment
    integer, intent(in), optional :: digits

    character(len=:), allocatable :: str
    integer :: i

    str = ''
    do i = 1, size(values)
       if (i > 1) str = str // ' '
       str = str // real_to_string(values(i), digits)
    enddo
    call write_key_value(key, str, comment)

  end subroutine log_list

!--------------------------------------------------------------------------------

  function real_to_string(x, digits, style) result(str)

    ! Write a real number compactly:
    !  - at most 'digits' significant digits (default 12), which hides floating-point noise
    !    such as 0.30000000000000004 but keeps values like 10.0000000001;
    !  - trailing zeros are dropped, but at least one digit follows the decimal point
    !    (910.0, 0.05, 2.4131e+04);
    !  - style = 'auto' (default): fixed notation for 1.e-6 <= |x| < 1.e6, scientific otherwise;
    !    style = 'fixed' or 'sci' forces one notation;
    !  - zero is written as 0.0.

    real(dp), intent(in) :: x
    integer, intent(in), optional :: digits
    character(len=*), intent(in), optional :: style
    character(len=:), allocatable :: str

    integer :: nd, ndec, e10, iexp, ipos
    real(dp) :: xr
    character(len=8) :: st
    character(len=64) :: buf, fmt
    character(len=:), allocatable :: mant
    logical :: use_sci

    nd = default_digits
    if (present(digits)) nd = max(1, min(digits, 17))
    st = 'auto'
    if (present(style)) st = style

    if (x == 0.0d0) then
       str = '0.0'
       return
    endif

    if (trim(st) == 'sci') then
       use_sci = .true.
    elseif (trim(st) == 'fixed') then
       use_sci = .false.
    else
       ! choose the notation from the value rounded to nd digits (e.g. 999999.9999999 rounds to 1.0e+06)
       write(fmt, '(a,i0,a,i0,a)') '(es', nd + 12, '.', nd - 1, 'e3)'
       write(buf, fmt) x
       read(buf, *) xr
       use_sci = (abs(xr) < fixed_min .or. abs(xr) >= fixed_max)
    endif

    if (.not. use_sci) then
       e10 = floor(log10(abs(x)))                ! power of 10 of the leading digit
       ndec = min(max(1, nd - (e10 + 1)), 40)    ! digits after the decimal point
       write(fmt, '(a,i0,a)') '(f60.', ndec, ')'
       write(buf, fmt) x
       str = trim(adjustl(buf))
       ! some compilers omit the zero before the decimal point
       if (str(1:1) == '.') str = '0' // str
       if (len(str) > 1) then
          if (str(1:2) == '-.') str = '-0' // str(2:)
       endif
       str = strip_zeros(str)
    else
       write(fmt, '(a,i0,a,i0,a)') '(es', nd + 12, '.', nd - 1, 'e3)'
       write(buf, fmt) x
       buf = adjustl(buf)
       ipos = index(buf, 'E')
       mant = strip_zeros(buf(1:ipos-1))
       read(buf(ipos+1:), *) iexp
       ! exponent with a sign and at least two digits (e-08, e+10, e-300)
       if (iexp < 0) then
          write(buf, '(a,"e-",i2.2)') mant, abs(iexp)
       else
          write(buf, '(a,"e+",i2.2)') mant, abs(iexp)
       endif
       if (abs(iexp) >= 100) then
          if (iexp < 0) then
             write(buf, '(a,"e-",i0)') mant, abs(iexp)
          else
             write(buf, '(a,"e+",i0)') mant, abs(iexp)
          endif
       endif
       str = trim(buf)
    endif

  end function real_to_string

!--------------------------------------------------------------------------------
! private procedures
!--------------------------------------------------------------------------------

  subroutine write_key_value(key, value, comment)

    ! Write '    key = value' with the key padded to key_width, and '  # comment' if present

    character(len=*), intent(in) :: key, value
    character(len=*), intent(in), optional :: comment

    character(len=:), allocatable :: line
    integer :: npad

    npad = max(0, key_width - len_trim(key))
    line = repeat(' ', indent_key) // trim(key) // repeat(' ', npad) // ' = ' // trim(value)
    if (present(comment)) then
       if (len_trim(comment) > 0) then
          npad = max(1, value_width - len_trim(value))
          line = line // repeat(' ', npad) // ' # ' // trim(comment)
       endif
    endif
    call write_log(line)

  end subroutine write_key_value

!--------------------------------------------------------------------------------

  function int_to_string(i) result(str)

    integer, intent(in) :: i
    character(len=:), allocatable :: str
    character(len=24) :: buf

    write(buf, '(i0)') i
    str = trim(buf)

  end function int_to_string

!--------------------------------------------------------------------------------

  function strip_zeros(s) result(str)

    ! Remove trailing zeros after a decimal point, keeping at least one digit after it

    character(len=*), intent(in) :: s
    character(len=:), allocatable :: str
    integer :: n, idot

    str = trim(s)
    idot = index(str, '.')
    if (idot == 0) return
    n = len(str)
    do while (n > idot + 1)
       if (str(n:n) /= '0') exit
       n = n - 1
    enddo
    str = str(1:n)

  end function strip_zeros

!--------------------------------------------------------------------------------

end module glissade_log_config
