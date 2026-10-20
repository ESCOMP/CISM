!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_lateral_melt.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glissade_lateral_melt

  use glide_types
  use glimmer_global, only: dp
  use glimmer_paramets, only: iulog, eps11
  use glimmer_log
  use glimmer_utils, only: point_diag

  use cism_parallel, only: this_rank, main_task, nhalo

  implicit none

  private
  public :: glissade_lateral_melt_solve, verbose_latmelt

  logical :: verbose_latmelt = .false.

contains

!-------------------------------------------------------------------------------

  subroutine glissade_lateral_melt_solve(model, lateral_melt)

    ! Compute lateral melting at the melt front, assuming either a constant melt rate
    !  or the parameterization of Rignot et al. (2016, GRL, doi:10.1002/2016GL068784).
    ! In principle, the melt front could include both grounded and floating ice,
    !  but runs to date (as of Oct. 2026) assume grounded ice only, with floating ice
    !  removed by the no-float calving scheme.

    use glimmer_physcon, only: scyr, scday
    use cism_parallel, only: parallel_type, parallel_halo, parallel_is_zero
    use glissade_masks, only: glissade_get_masks, glissade_melt_front_mask
    use glissade_calving, only: compute_calving_front_length, calve_ice_basic

    implicit none

    ! input/output arguments

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! Passing in lateral melt in addition to model saves some typing below
    type(glide_lateral_melt), intent(inout) :: lateral_melt

    ! --- Local variables ---

    integer, dimension(model%general%ewn, model%general%nsn) :: &
         ice_mask,                & ! = 1 if ice is present
         floating_mask,           & ! = 1 if ice is present and floating
         land_mask,               & ! = 1 if topg - eus >= 0
         ocean_mask                 ! = 1 if ice is absent and topg - eus < 0

    real(dp), dimension(model%general%ewn, model%general%nsn) :: &
         thck_submerged,          & ! thickness of submerged ice (m)
         mf_length,               & ! melt front length for each cell (m)
         surface_runoff             ! surface meltwater runoff (kg/m^2/s)

    integer :: nx, ny               ! horizontal grid dimensions
    real(dp) :: dx, dy              ! cell dimensions in x and y directions (m)
    real(dp) :: dt                  ! timestep (s)
    integer :: itest, jtest, rtest  ! coordinates of diagnostic point
    integer :: i, j

    type(parallel_type) :: parallel   ! info for parallel communication

    nx = model%general%ewn
    ny = model%general%nsn
    dx = model%numerics%dew
    dy = model%numerics%dns
    dt = model%numerics%dt

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    parallel = model%parallel

    call parallel_halo(model%geometry%thck, parallel)

    ! Compute some masks

    call glissade_get_masks(&
         nx,                ny,            &
         parallel,                         &
         model%geometry%thck,              &  ! m
         model%geometry%topg,              &  ! m
         model%climate%eus,                &  ! m
         eps11,                            &  ! thklim (m) = eps11
         ice_mask,                         &
         floating_mask = floating_mask,    &
         ocean_mask = ocean_mask,          &
         land_mask = land_mask)

    call glissade_melt_front_mask(&
         nx,                ny,            &
         ice_mask,                         &
         land_mask,         ocean_mask,    &
         lateral_melt%melt_front_mask)

    call parallel_halo(lateral_melt%melt_front_mask, parallel)

    ! Compute the melt front length for each cell based on its number of ocean neighbors
    ! (i.e., enhanced lateral melting for cells with 2 or 3 ocean neighbors).
    ! Call the subroutine used to compute the calving front length based on calving_front_mask.

    call compute_calving_front_length(&
         nx,                ny,            &
         dx,                dy,            &
         itest,   jtest,    rtest,         &
         lateral_melt%melt_front_mask,     &
         ocean_mask,                       &
         mf_length)

    ! Compute the submerged ice thickness.
    ! It is set to the negative of the topography for marine-grounded ice,
    !  and zero for land-grounded ice.

    thck_submerged = model%calving%thck_effective*(rhoi/rhoo)
    thck_submerged = min(thck_submerged, max(model%climate%eus - model%geometry%topg, 0.0d0))

    if (verbose_latmelt) then
       call point_diag(model%geometry%thck, 'thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(model%geometry%topg, 'topg (m)', itest, jtest, rtest, 7, 7)
       call point_diag(thck_submerged, 'thck_submerged (m)', itest, jtest, rtest, 7, 7)
       call point_diag(lateral_melt%melt_front_mask, 'melt_front_mask', itest, jtest, rtest, 7, 7)
       call point_diag(ocean_mask, 'ocean_mask', itest, jtest, rtest, 7, 7)
       call point_diag(mf_length, 'mf_length (m)', itest, jtest, rtest, 7, 7)
    endif

    if (model%options%which_lateral_melt == LATERAL_MELT_CONSTANT) then

       call constant_lateral_melt(&
            nx,                   ny,              &
            dx,                   dy,              &
            dt,                                    &  ! s
            itest,    jtest,      rtest,           &
            lateral_melt%melt_front_mask,          &
            lateral_melt%melt_rate_const/scyr,     &  ! m/s
            thck_submerged,                        &  ! m
            mf_length,                             &  ! m
            model%geometry%thck,                   &  ! m
            lateral_melt%latmelt_thck)                ! m

       if (verbose_latmelt .and. main_task) then
          write(iulog,*) 'Called constant_lateral_melt, melt_rate_const (m/yr) =', &
               lateral_melt%melt_rate_const
       endif

    elseif (model%options%which_lateral_melt == LATERAL_MELT_ISMIP) then

       if (lateral_melt%runoff_from_ablation) then
          ! compute from the negative portion of acab
          surface_runoff = max(-1.0d0*rhoi*model%climate%acab_applied, 0.0d0)   ! convert m/s ice to kg/m^2/s w.e.

          call surface_runoff_to_discharge(&
               nx,                  ny,             &
               dx,                  dy,             &  ! m
               parallel,                            &
               model%ocean_data%nbasin,             &
               model%ocean_data%basin_number,       &
               ice_mask,                            &
               surface_runoff,                      &  ! kg/m^2/s
               thck_submerged,                      &  ! m
               mf_length,                           &  ! m
               lateral_melt%runoff_applied)            ! kg/m^2/s

       endif

       if (lateral_melt%thermal_forcing_avg_3d_to_2d) then

          if (parallel_is_zero(model%ocean_data%thermal_forcing)) then
             call write_log('Missing 3D thermal forcing for lateral melt', GM_FATAL)
          endif

          call thermal_forcing_avg_3d_to_2d(&
               nx,                  ny,             &
               model%ocean_data%nzocn,              &
               model%ocean_data%zocn,               &
               model%ocean_data%thermal_forcing,    &  ! K
               lateral_melt%ztop_tfavg,             &  ! m
               lateral_melt%zbot_tfavg,             &  ! m
               lateral_melt%thermal_forcing_applied)   ! K

       else   ! should have received the 2D thermal forcing directly (from input file or coupler)

          if (parallel_is_zero(model%lateral_melt%thermal_forcing_applied)) then
             call write_log('Missing thermal_forcing_applied for lateral melt', GM_FATAL)
          endif

       endif

       call ismip_lateral_melt(&
            nx,                 ny,                &
            dx,                 dy,                &
            dt,                                    &  ! s
            itest,   jtest,     rtest,             &
            lateral_melt%melt_front_mask,          &
            lateral_melt%melt_factor,              &
            lateral_melt%runoff_applied,           &  ! kg/m2/s
            lateral_melt%thermal_forcing_applied,  &  ! K
            thck_submerged,                        &  ! m
            mf_length,                             &  ! m
            model%geometry%thck,                   &  ! m
            lateral_melt%latmelt_thck)                ! m

       if (verbose_latmelt) then
          if (main_task) write(iulog,*) 'Called ismip_lateral_melt, melt_factor =', lateral_melt%melt_factor
          call point_diag(lateral_melt%runoff_applied*scday, 'runoff_applied (kg/m^2/day)', &
               itest, jtest, rtest, 7, 7)
          call point_diag(lateral_melt%thermal_forcing_applied, 'thermal_forcing_applied (deg C)', &
               itest, jtest, rtest, 7, 7)
       endif

    endif  ! which_lateral_melt

    if (verbose_latmelt) then
       call point_diag(lateral_melt%latmelt_thck, 'latmelt_thck (m)', itest, jtest, rtest, 7, 7)
    endif

    call parallel_halo(model%geometry%thck, parallel)
    call parallel_halo(lateral_melt%latmelt_thck, parallel)

    ! Melting may have brought some cells to flotation
    ! In the case CALVING_FLOAT_ZERO (which is the usual case when using a lateral melt parameterization),
    !  remove floating ice so that none exists at the end of the timestep.

    if (model%options%whichcalving == CALVING_FLOAT_ZERO) then

       call calve_ice_basic(&
            nx,             ny,             &
            model%options%whichcalving,     &
            model%options%calving_domain,   &
            parallel,                       &
            model%calving,                  &  ! calving derived type
            itest,  jtest,  rtest,          &
            dt,                             &  ! s
            eps11,                          &  ! thcklim (m) = eps11
            model%geometry%thck,            &  ! m
            model%isostasy%relx,            &  ! m
            model%geometry%topg,            &  ! m
            model%climate%eus)                 ! m

    endif  ! CALVING_FLOAT_ZERO

    if (verbose_latmelt) then
       call point_diag(model%geometry%thck, 'New thck (m)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine glissade_lateral_melt_solve

!-------------------------------------------------------------------------------

  subroutine constant_lateral_melt(&
       nx,                 ny,           &
       dx,                 dy,           &
       dt,                               &  ! s
       itest,   jtest,     rtest,        &
       melt_front_mask,                  &
       melt_rate_const,                  &  ! m/s
       thck_submerged,                   &  ! m
       mf_length,                        &  ! m
       thck,                             &  ! m
       latmelt_thck)                        ! m

    ! Apply lateral melt horizontally based on a prescribed constant melt rate

    use glimmer_physcon, only: rhoi, rhoo, scyr

    ! input/output arguments
    
    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt                        ! time step (s)

    integer, dimension(nx,ny), intent(in)  ::  &
         melt_front_mask           ! = 1 where ice is marine-grounded or floating
                                   ! and borders at least one ocean cell, else = 0

    real(dp), intent(in) :: &
         melt_rate_const           ! prescribed constant melt rate (m/s)

    real(dp), dimension(nx,ny), intent(in) :: &
         thck_submerged,         & ! effective thickness (m) of submerged ice
         mf_length                 ! length of melt front in each grid cell (m)

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck                      ! ice thickness (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         latmelt_thck              ! thickness reduction (m) due to lateral melt

    ! local variables

    integer :: i, j

!!    real(dp) :: &
!!         m_sr                    ! horizontal melting rate in m/yr calculated from Slater ISMIP6 melt approach

    ! Initialize

    latmelt_thck = 0.0d0

    ! Loop over locally owned cells
    ! Melt occurs only in MF cells: marine ice-filled cells with one or more ocean neighbors.
    
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (melt_front_mask(i,j) == 1) then

             ! Apply an absolute “horizontal” melt rate - i.e., the melt rate perpendicular
             ! to the face of an approximately vertical calving front
!!             m_sr = frontal_melt_rate ! m/s
!!             latmelt_dthck(i,j) = min((m_sr*dt * thck_effective(i,j) * mf_length(i,j)) / (dx*dy), thck(i,j))

             !WHL - commented out m_sr; just using the prescribed melt rate with thck_submerged
             !TODO - Question: What happens to the ice above sea level? Should it collapse?
             latmelt_thck(i,j) = melt_rate_const * dt * thck_submerged(i,j) * mf_length(i,j) / (dx*dy)

             ! limit latmelt_thck and adjust thck
             latmelt_thck (i,j) = min(latmelt_thck(i,j), thck(i,j))
             thck(i,j) = thck(i,j) - latmelt_thck(i,j)

             if (verbose_latmelt) then
                if (this_rank == rtest .and. i == itest .and. j == jtest) then
                   write(iulog,*) 'Constant lateral melting, rank, i, j:', this_rank, i, j
                   write(iulog,*) 'melt rate (m/yr), Hsub, mf_length, latmelt_thck, new thck:', &
                        melt_rate_const*scyr, thck_submerged(i,j), mf_length(i,j), latmelt_thck(i,j), thck(i,j)
                endif
             endif

          endif   ! melt_front_mask = 1
       enddo   ! i
    enddo   ! j

  end subroutine constant_lateral_melt

!-------------------------------------------------------------------------------

  subroutine ismip_lateral_melt(&
       nx,                 ny,           &
       dx,                 dy,           &
       dt,                               &  ! s
       itest,   jtest,     rtest,        &
       melt_front_mask,                  &
       melt_factor,                      &
       subglacial_discharge,             &  ! m/s
       thermal_forcing_2d,               &  ! K
       thck_submerged,                   &  ! m
       mf_length,                        &  ! m
       thck,                             &  ! m
       melt_thck)                           ! m

    ! Apply lateral melt horizontally as a function of subglacial discharge and thermal forcing.
    ! Based on the parameterization of X.

    use glimmer_physcon, only: rhow, scday, scyr

    ! input/output arguments
    
    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt                        ! time step (s)

    !TODO - Sufficient to pass mf_length?
    integer, dimension(nx,ny), intent(in)  ::  &
         melt_front_mask           ! = 1 where floating or marine-grounded ice borders the ocean

    real(dp), intent(in) :: &
         melt_factor               ! multiplier for Rignot melt parameterisation (unitless)

    real(dp), dimension(nx,ny), intent(in) :: &
         subglacial_discharge,   & ! subglacial meltwater discharge (kg/m^2/s)
         thermal_forcing_2d,     & ! average thermal forcing over some depth range (K)
         thck_submerged,         & ! effective thickness (m) of submerged ice
         mf_length                 ! length of melt front in each grid cell (m)

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck                      ! ice thickness (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         melt_thck                 ! thickness reduction (m) due to lateral melt

    ! local variables

    integer :: i, j

    real(dp) :: &
         q_sr,                   & ! runoff in Rignot calculation in m/d
         tf_sr,                  & ! thermal forcing in deg C
         m_sr                      ! melting rate in m/yr calculated from Slater ISMIP6 melt approach

    ! Initialize

    melt_thck = 0.0d0

    ! Loop over locally owned cells
    ! Melt occurs only in MF cells: marine ice-filled cells with one or more ocean neighbors.
    
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (melt_front_mask(i,j) == 1) then   !TODO - replace with mf_length > 0.0

             ! Apply an absolute “horizontal” melt rate - i.e., the melt rate perpendicular
             !  to the face of an approximately vertical melt front
             ! Parameterize melt as a function of subglacial discharge, thermal forcing and effective thickness
             ! See Rignot et al. 2016 or ISMIP6 melt forcing approach.

             tf_sr = thermal_forcing_2d(i,j)  ! 2d thermal forcing [degC]
             q_sr = (subglacial_discharge(i,j)/rhow) * scday  ! subglacial_discharge converted from kg/m^2/s to m/day

             ! Rignot et al. 2016; formulated in m/day, convert to m/s and multiply by a prescribed melt_factor
             m_sr = melt_factor * (3.0d-4 * thck_submerged(i,j) * q_sr**0.39d0 + 0.15d0) * tf_sr**1.18d0 * 365.d0/scyr

             ! calculate applied thickness change
             melt_thck(i,j) = m_sr * dt * thck_submerged(i,j) * mf_length(i,j) / (dx*dy)

             ! limit melt_thck and adjust thck
             melt_thck (i,j) = min(melt_thck(i,j), thck(i,j))
             thck(i,j) = thck(i,j) - melt_thck(i,j)

             if (verbose_latmelt) then
                if (this_rank == rtest .and. i == itest .and. j == jtest) then
                   write(iulog,*) 'ISMIP6 lateral melting, rank, i, j:', this_rank, i, j
                   write(iulog,*) 'Hsub, m_sr (m/yr), mf_length, melt_thck:', &
                        thck_submerged(i,j), m_sr*365.d0, mf_length(i,j), melt_thck(i,j)
                endif
             endif

          endif
       enddo
    enddo

  end subroutine ismip_lateral_melt

!-------------------------------------------------------------------------------

  subroutine surface_runoff_to_discharge(&
       nx,             ny,                &
       dx,             dy,                &
       parallel,                          &
       nbasin,         basin_number,      &
       ice_mask,                          &
       surface_runoff,                    &
       thck_submerged,                    &
       mf_length,                         &
       runoff_applied)

    ! Given the surface runoff, compute basin-scale subglacial discharge ('runoff_applied').
    ! This as an input to the lateral melt parameterization.

    use glimmer_physcon, only: rhoi, rhow, rhoo
    use glissade_utils, only: glissade_basin_sum

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny                        !> number of grid cells in each dimension

    real(dp), intent(in) :: &
         dx, dy                        !> grid cell size (m)

    type(parallel_type), intent(in) :: &
         parallel                      !> info for parallel communication

    integer, intent(in) :: nbasin      !> number of basins

    integer, dimension(nx,ny), intent(in) :: &
         basin_number,               & !> basin number for each grid cell
         ice_mask                      !> = 1 where ice is present, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         surface_runoff,             & !> surface runoff (kg/m^2/s)
         thck_submerged,             & !> effective thickness (m) of submerged ice
         mf_length                     ! length of melt front in each grid cell (m)

    ! discharge per unit area at the submerged ice front; nonzero only for cells at the front
    real(dp), dimension(nx,ny), intent(out) :: &
         runoff_applied                !> subglacial discharge applied to lateral melting (kg/m2/s)

    ! local variables

    integer :: i, j, nb

    real(dp), dimension(nx,ny) :: &
         rmask                         ! real mask for basin sums

    real(dp), dimension(nbasin) :: &
         runoff_sum_basin,           & ! runoff summed over each basin (m^3/s)
         area_submerged_sum_basin,   & ! submerged area summed over each basin (m^2)
         runoff_applied_basin          ! discharge per unit front area in each basin (kg/m^2/s)

    ! Sum the runoff over each basin

    rmask = 1.0d0 * ice_mask

    call glissade_basin_sum(&
         nx,         ny,                &
         parallel,                      &
         nbasin,     basin_number,      &
         rmask,                         &
         surface_runoff*dx*dy,          & ! kg/s
         runoff_sum_basin)                ! kg/s

    ! Sum the submerged ice area over each basin

    call glissade_basin_sum(&
         nx,         ny,                &
         parallel,                      &
         nbasin,     basin_number,      &
         rmask,                         &
         mf_length*thck_submerged,      & ! m^2
         area_submerged_sum_basin)        ! m^2

    ! Compute the average subglacial discharge in each basin as the runoff per unit submerged melt front area

    runoff_applied_basin = 0.0d0
    do nb = 1, nbasin
       if (area_submerged_sum_basin(nb) > 0.0d0) then
          ! Divide total basin runoff (kg/s) by total basin submerged area (m^2)
          runoff_applied_basin(nb) = runoff_sum_basin(nb) / area_submerged_sum_basin(nb)  ! kg/m^2/s
       endif
    enddo

    ! Copy the basin-average subglacial discharge to cells at the melt front

    runoff_applied = 0.0d0
    do j = 1, ny
       do i = 1, nx
          nb = basin_number(i,j)
          if (nb >= 1 .and. nb <= nbasin .and. mf_length(i,j) > 0.0d0) then
             runoff_applied(i,j) = runoff_applied_basin(nb)
          endif
       enddo   ! i
    enddo   ! j

  end subroutine surface_runoff_to_discharge

!-------------------------------------------------------------------------------

  subroutine thermal_forcing_avg_3d_to_2d(&
       nx,              ny,      &
       nzocn,                    &
       zocn,                     &
       thermal_forcing,          &
       ztop,            zbot,    &
       thermal_forcing_2d)

    ! Average the 3D thermal forcing over a prescribed depth range.
    ! The default range is -200 m to -500 m.
    ! Note: This could be converted to a utility model with an arbitrary 3d field as input.
    !       In the lateral melt module for now since not yet used for fields other than TF.

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny                    !> number of grid cells in each dimension

    integer, intent(in) :: &
         nzocn                     !> number of ocean levels

    real(dp), dimension(nzocn), intent(in) :: &
         zocn                      !> ocean levels (m) where forcing is provided, negative below sea level

    real(dp), dimension(nzocn,nx,ny), intent(in) :: &
         thermal_forcing           !> thermal forcing field at ocean levels

    real(dp), intent(in)  :: &
         ztop, zbot                !> top and bottom of depth range (m), negative below sea level
                                   !> default values are -200 m and -500 m

    real(dp), dimension(nx,ny), intent(out) :: &
         thermal_forcing_2d        !> average thermal forcing over some depth range

    ! local variables

    integer :: i, j, k

    real(dp) :: &
         layer_frac,             & ! fraction of layer within the depth range
         dlayer,                 & ! layer thickness
         thermal_forcing_layer     ! thermal forcing in the layer

    real(dp), dimension(0:nzocn) :: zbnd   ! depths of layer boundaries

    if (ztop >= 0.0d0 .or. zbot >= 0.0d0) then
       call write_log('Error, average_thermal_forcing, ztop and zbot must be < 0', GM_FATAL)
    endif

    ! Estimate the boundaries between ocean layers
    ! Note: k = 1 is the top level, and zocn becomes more negative with increasing k.
    ! zocn(k) is the depth and the middle of layer k, and zbnd(k) is the depth at the bottom of layer k.
    ! For uniform layers, the spacing between boundaries is the same as the spacing between layers.

    zbnd(0) = 0.0d0
    do k = 1, nzocn-1
       zbnd(k) = 0.5d0 * (zocn(k) + zocn(k+1))
    enddo
    zbnd(nzocn) = zocn(k) - 0.5d0*(zocn(nzocn-1) - zocn(nzocn))

    if (verbose_latmelt .and. main_task) then
       write(iulog,*) 'ocean layers, k, zocn, zbnd:'
       do k = 1, nzocn
          write(iulog,*) k, zocn(k), zbnd(k)
       enddo
    endif

    ! Average the thermal forcing over the specified depth range

    thermal_forcing_2d = 0.0d0

    do k = 1, nzocn
       if (zbnd(k) < ztop .and. zbnd(k-1) > zbot) then   ! include this layer in the average
          if (zbnd(k-1) > ztop) then  ! part of the layer is above the range
             layer_frac = (ztop - zbnd(k)) / (zbnd(k-1) - zbnd(k))
          elseif (zbnd(k) < zbot) then  ! part of the layer is below the range
             layer_frac = (zbnd(k-1) - zbot) / (zbnd(k-1) - zbnd(k))
          else
             layer_frac = 1.0d0  ! the entire layer is within the range
          endif
          dlayer = zbnd(k-1) - zbnd(k)
          if (main_task .and. verbose_latmelt) write(iulog,*) 'k, dlayer, layer_frac:', k, dlayer, layer_frac

          do j = 1, ny
             do i = 1, nx
                ! Accumulate the forcing for this layer
                thermal_forcing_2d(i,j) = thermal_forcing_2d(i,j) + thermal_forcing_layer*dlayer*layer_frac
             enddo
          enddo
       endif
    enddo

    ! Divide by the depth range
    where (thermal_forcing_2d > 0.0d0)
       thermal_forcing_2d = thermal_forcing_2d / (ztop - zbot)
    endwhere

  end subroutine thermal_forcing_avg_3d_to_2d

!-------------------------------------------------------------------------------

end module glissade_lateral_melt

!-------------------------------------------------------------------------------
