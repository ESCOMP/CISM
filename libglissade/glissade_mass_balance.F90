!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_mass_balance.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains drivers for handling the surface and basal mass balance.
!
! Author: William Lipscomb
!         NSF National Center for Atmospheric Research
!         Climate and Global Dynamics Laboratory
!         Boulder, CO 80303
!         USA
!         <lipscomb@ucar.edu>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
  module glissade_mass_balance

    use glimmer_global, only: dp
    use glimmer_log
    use glide_types
    use cism_parallel, only: this_rank, main_task, nhalo, &
         parallel_reduce_max

    implicit none
    save
    private

    public :: glissade_mass_balance_init, glissade_mass_balance_solve
!              glissade_overwrite_acab_mask, glissade_overwrite_acab,  &
!              glissade_add_2d_anomaly, glissade_add_3d_anomaly

!    logical, parameter ::     &
!         conservation_check = .true. ! if true, check global conservation

!=======================================================================

  contains

!=======================================================================

  subroutine glissade_mass_balance_init(model)

    ! Initialize some fields related to the surface mass balance

    use glimmer_paramets, only: eps11
    use glimmer_physcon, only: rhow, rhoi
    use glimmer_scales, only: scale_acab

    ! input/output arguments

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! local variables

    real(dp) :: local_maxval, global_maxval
    character(len=100) :: message
      
    ! Initialize acab, if SMB (with different units) was read in
    if (model%options%smb_input == SMB_INPUT_MMYR_WE) then
       ! Convert units from mm/yr w.e. to m/yr ice, then convert to model units
       model%climate%acab(:,:) = (model%climate%smb(:,:) * (rhow/rhoi)/1000.d0) / scale_acab
       !WHL - debug
       if (main_task) write(6,*) 'Setting acab, m/yr ice'
    endif

    ! Initialize artm_corrected.  This is equal to artm, plus any prescribed temperature anomaly.
    model%climate%artm_corrected(:,:) = model%climate%artm(:,:)
    
    if (model%options%enable_artm_anomaly) then
       ! Check whether artm_anomaly was read from an external file.
       ! If so, then use this field as the anomaly.
       ! If not, then set artm_anomaly = artm_anomaly_constant everywhere.
       ! Note: The artm_anomaly field does not change during the run,
       !       but it is possible to ramp up the anomaly using artm_anomaly_timescale.

       local_maxval = maxval(abs(model%climate%artm_anomaly))
       global_maxval = parallel_reduce_max(local_maxval)
       if (global_maxval < eps11) then
          model%climate%artm_anomaly = model%climate%artm_anomaly_const
          write(message,*) &
               'Setting artm_anomaly = constant value (degC):', model%climate%artm_anomaly_const
          call write_log(trim(message))
       else
          if (model%options%is_restart == NO_RESTART) then
             call write_log('Setting artm_anomaly from external file')
          endif
       endif
    endif
    !TODO - Write a short utility function to compute global_maxval of any field.
    !TODO - Repeat for snow and precip anomalies

    ! If acab is to be overwritten for some cells, then set overwrite_acab_mask = 1 for these cells.
    ! We can overwrite the input acab with a fixed value (typically negative) where
    ! (1) the input acab = 0 at initialization, or
    ! (2) the input thck <= overwrite_acab_minthck at initialization
    ! Note: This option is designed for standalone runs, and should be used with caution for coupled runs.
    !       On restart, overwrite_acab_mask is read from the restart file.
    
    if (model%climate%overwrite_acab_value /= 0 .and. model%options%is_restart == NO_RESTART) then

       call set_overwrite_acab_mask(&
            model%options%overwrite_acab,          &
            model%climate%acab,                    &
            model%geometry%thck,                   &
            model%climate%overwrite_acab_minthck,  &
            model%climate%overwrite_acab_mask)

    endif

  end subroutine glissade_mass_balance_init

!=======================================================================

  subroutine glissade_mass_balance_solve

  end subroutine glissade_mass_balance_solve

!=======================================================================

  subroutine set_overwrite_acab_mask(&
       overwrite_acab,         &
       acab,                   &
       thck,                   &
       overwrite_acab_minthck, &
       overwrite_acab_mask)

    use glide_types

    ! If overwrite_acab /=0 , then set overwrite_acab_mask = 1 for grid cells
    !  where acab is to be overwritten.  Currently, three options are supported:
    ! (1) Overwrite acab where the input acab = 0 at initialization
    ! (2) Overwrite acab where the input thck <= overwrite_acab_minthck at initialization
    ! (3) Overwrite acab based on an input mask
    !
    ! Note: This subroutine should be called only on initialization, not on restart.

    integer, intent(in) ::  &
         overwrite_acab           !> option for overwriting acab

    real(dp), dimension(:,:), intent(in) ::  &
         acab,                  & !> ice surface mass balance (model units)
         thck                     !> ice thickness (model units)

    real(dp), intent(in) ::  &
         overwrite_acab_minthck   !> overwrite acab where thck <= overwrite_acab_minthck (model units)

    integer, dimension(:,:), intent(out) ::  &
         overwrite_acab_mask      !> = 1 where acab is overwritten, else = 0

    integer :: ewn, nsn
    integer :: i, j
    integer :: max_mask_local, max_mask_global

    ewn = size(overwrite_acab_mask,1)
    nsn = size(overwrite_acab_mask,2)

    if (overwrite_acab == OVERWRITE_ACAB_ZERO_ACAB) then

       do j = 1, nsn
          do i = 1, ewn

             if (acab(i,j) == 0.0d0) then
                overwrite_acab_mask(i,j) = 1
             else
                overwrite_acab_mask(i,j) = 0
             endif

          enddo
       enddo

    elseif (overwrite_acab == OVERWRITE_ACAB_THCKMIN) then

       do j = 1, nsn
          do i = 1, ewn

             ! Note the '<='.  If overwrite_acab_minthck = 0.d0, only ice-free cells are overwritten.
             if (thck(i,j) <= overwrite_acab_minthck) then
                overwrite_acab_mask(i,j) = 1
             else
                overwrite_acab_mask(i,j) = 0
             endif

          enddo
       enddo

    elseif (overwrite_acab == OVERWRITE_ACAB_INPUT_MASK) then

       ! Make sure a mask was read in with some nonzero values
       ! If not, then write a warning

       max_mask_local = maxval(overwrite_acab_mask)
       max_mask_global = parallel_reduce_max(max_mask_local)
       if (main_task) then
          print*, 'rank, max_mask_local, max_mask_global:', &
               this_rank, max_mask_local, max_mask_global
       endif
       if (max_mask_global == 1) then
          ! continue
       elseif (max_mask_global == 0) then
          call write_log('Using overwrite_acab_mask without any values > 0', GM_WARNING)
       else
          call write_log('Using overwrite_acab_mask with values other than 0 and 1', GM_FATAL)
       endif

    endif  ! overwrite_acab

  end subroutine set_overwrite_acab_mask

!=======================================================================

  end module glissade_mass_balance

!=======================================================================
