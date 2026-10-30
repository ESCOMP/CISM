!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   glissade_types.F90 - part of the Community Ice Sheet Model (CISM)
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

module glissade_types

  ! Derived types, option constants and parameters for the Glissade dycore.
  !
  ! For now, this module re-exports everything in glide_types, so that
  !  libglissade modules can use glissade_types instead of glide_types.
  ! It also makes glide_global_type available as glissade_global_type, the
  !  name used by Glissade-specific code generated from ncdf_template.F90.in.
  ! The Glissade types will be defined here when Glissade gets its own
  !  model instance type.

  use glide_types
  use glide_types, only: glide_global_type, glissade_global_type => glide_global_type

  implicit none

end module glissade_types
