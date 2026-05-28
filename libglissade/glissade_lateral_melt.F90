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
  public :: glissade_lateral_melt_constant, glissade_lateral_melt_ismip, &
       glissade_thermal_forcing_avg_3d_to_2d, glissade_subglacial_discharge

  public :: verbose_latmelt

  logical, parameter :: verbose_latmelt = .true.

contains

!-------------------------------------------------------------------------------

  subroutine glissade_lateral_melt_constant(&
       nx,                 ny,           &
       dx,                 dy,           &
       dt,                 time,         &  ! s
       itest,   jtest,     rtest,        &
       melt_front_mask,                  &
       melt_rate_const,                  &  ! m/s
       thck_submerged,                   &  ! m
       mf_length,                        &  ! m
       latmelt_dthck)                       ! m

    ! Apply lateral melt horizontally based on a prescribed constant melt rate

    use glimmer_physcon, only: rhoi, rhoo, scyr

    ! input/output arguments
    
    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt,                     & ! time step (s)
         time                      ! elapsed time (s) of model run

    integer, dimension(nx,ny), intent(in)  ::  &
         melt_front_mask           ! = 1 where ice is marine-grounded or floating
                                   ! and borders at least one ocean cell, else = 0

    real(dp), intent(in) :: &
         melt_rate_const           ! prescribed constant melt rate (m/yr)

    real(dp), dimension(nx,ny), intent(in) :: & 
         thck_submerged,         & ! effective thickness (m) of submerged ice
         mf_length                 ! length of melt front in each grid cell (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         latmelt_dthck             ! thickness reduction (m) due to lateral melt

    ! local variables

    integer :: i, j

!!    real(dp) :: &
!!         m_sr                    ! horizontal melting rate in m/yr calculated from Slater ISMIP6 melt approach

    ! Initialize

    latmelt_dthck = 0.0d0

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
             latmelt_dthck(i,j) = melt_rate_const * dt * thck_submerged(i,j) * mf_length(i,j) / (dx*dy)

             if (verbose_latmelt) then
                if (this_rank == rtest .and. i == itest .and. j == jtest) then
                   write(iulog,*) 'Constant lateral melting, rank, i, j:', this_rank, i, j
                   write(iulog,*) 'melt rate (m/yr), Hsub, mf_length, latmelt_dthck:', &
                        melt_rate_const*scyr, thck_submerged(i,j), mf_length(i,j), latmelt_dthck(i,j)
                endif
                call point_diag(latmelt_dthck, 'lateral melt dthck', itest, jtest, rtest, 7, 7)
             endif

          endif   ! melt_front_mask = 1
       enddo   ! i
    enddo   ! j

  end subroutine glissade_lateral_melt_constant

!-------------------------------------------------------------------------------

  subroutine glissade_lateral_melt_ismip(&
       nx,                 ny,           &
       dx,                 dy,           &
       dt,                 time,         &  ! s
       itest,   jtest,     rtest,        &
       melt_front_mask,                  &
       melt_factor,                      &
       subglacial_discharge,             &  ! m/s
       thermal_forcing_2d,               &  ! K
       thck_submerged,                   &  ! m
       mf_length,                        &  ! m
       latmelt_dthck)                       ! m

    ! Apply lateral melt horizontally as a function of subglacial discharge and thermal forcing.
    ! Based on the parameterization of X.

    use glimmer_physcon, only: rhoi, rhoo, scday, scyr

    ! input/output arguments
    
    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt,                     & ! time step (s)
         time                      ! elapsed time (s) of model run

    integer, dimension(nx,ny), intent(in)  ::  &
         melt_front_mask           ! = 1 where floating or marine-grounded ice borders the ocean

    real(dp), intent(in) :: &
         melt_factor               ! multiplier for Rignot melt parameterisation

    real(dp), dimension(nx,ny), intent(in) :: & 
         subglacial_discharge,   & ! subglacial meltwater discharge (m/s)
         thermal_forcing_2d,     & ! average thermal forcing over some depth range (K)
         thck_submerged,         & ! effective thickness (m) of submerged ice
         mf_length                 ! length of melt front in each grid cell (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         latmelt_dthck             ! thickness reduction (m) due to lateral melt

    ! local variables

    integer :: i, j

    real(dp) :: &
         q_sr,                   & ! runoff in Rignot calculation in m/d
         tf_sr,                  & ! thermal forcing in deg C
         m_sr                      ! melting rate in m/yr calculated from Slater ISMIP6 melt approach

    ! Initialize

    latmelt_dthck = 0.0d0

    ! Loop over locally owned cells
    ! Melt occurs only in MF cells: marine ice-filled cells with one or more ocean neighbors.
    
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (melt_front_mask(i,j) == 1) then

             ! Apply an absolute “horizontal” melt rate - i.e., the melt rate perpendicular
             !  to the face of an approximately vertical calving front
             ! Parameterize melt as a function of subglacial discharge, thermal forcing and effective thickness
             ! See Rignot et al. 2016 or ISMIP6 melt forcing approach.

             tf_sr = thermal_forcing_2d(i,j) ! 2d thermal forcing [degC]
             q_sr = subglacial_discharge(i,j) * scday ! subglacial_discharge passed in m/s; for Rignot equation convert to [m/d]

             ! Rignot et al. 2016; formulated in m/d, converted to m/s. Mulitplier frontal_melt_factor as proposed for ISMIP7
             m_sr = melt_factor * (3.0d-4 * thck_submerged(i,j) * q_sr**0.39d0 + 0.15d0) * tf_sr**1.18d0 * 365.d0/scyr

             ! calculate applied thickness change
             latmelt_dthck(i,j) = m_sr * dt * thck_submerged(i,j) * mf_length(i,j) / (dx*dy)

             if (verbose_latmelt) then
                if (this_rank == rtest .and. i == itest .and. j == jtest) then
                   write(iulog,*) 'ISMIP6 lateral melting, rank, i, j:', this_rank, i, j
                   write(iulog,*) 'Hsub, mf_length, latmelt_dthck:', &
                        thck_submerged(i,j), mf_length(i,j), latmelt_dthck(i,j)
                endif
                call point_diag(latmelt_dthck, 'lateral melt dthck', itest, jtest, rtest, 7, 7)
             endif

          endif
       enddo
    enddo

  end subroutine glissade_lateral_melt_ismip

!-------------------------------------------------------------------------------

  subroutine glissade_subglacial_discharge(&
       nx,             ny,                &
       dx,             dy,                &
       parallel,                          &
       nbasin,         basin_number,      &
       ice_mask,                          &
       acab,                              &
       thck_submerged,                    &
       mf_length,                         &
       subglacial_discharge)

    ! Compute basin-scale subglacial discharge.
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
         acab,                       & !> applied accumulation/ablation (m/s)
         thck_submerged,             & !> effective thickness (m) of submerged ice
         mf_length                     ! length of melt front in each grid cell (m)

    ! Note: To match the ISMIP6 input file units (kg/m^2/s), multiply by rhow = 1000 kg/m^3
    real(dp), dimension(nx,ny), intent(out) :: &
         subglacial_discharge          !> subglacial meltwater discharge for lateral melting (m/s)

    ! local variables

    integer :: i, j, nb

    real(dp), dimension(nx,ny) :: &
         runoff,                     & ! subglacial runoff, estimated as a function of surface ablation (m^3/s)
         rmask                         ! real mask for basin sums

    real(dp), dimension(nbasin) :: &
         runoff_sum_basin,           & ! runoff summed over each basin (m^3/s)
         area_submerged_sum_basin      ! submerged area summed over each basin (m^2)

    ! Estimate the surface runoff reaching the bed of the ice sheet.
    ! Assume that all the surface ablation (acab, wherever acab < 0) reaches the bed.

    runoff = max(-1.0d0*acab*(rhoi/rhow), 0.0d0) * (dx*dy)   ! m^3/s

    ! Sum the runoff over each basin

    rmask = 1.0d0 * ice_mask

    call glissade_basin_sum(&
         nx,         ny,                &
         parallel,                      &
         nbasin,     basin_number,      &
         rmask,                         &
         runoff,                        & ! m^3/s
         runoff_sum_basin)                ! m^3/s

    ! Sum the submerged ice area over each basin

    call glissade_basin_sum(&
         nx,         ny,                &
         parallel,                      &
         nbasin,     basin_number,      &
         rmask,                         &
         mf_length*thck_submerged,      & ! m^2
         area_submerged_sum_basin)        ! m^2

    ! Compute the subglacial discharge as the runoff per unit submerged melt front area.
    ! This discharge is uniform across each basin.
    ! The result is similar to what ISMIP6 provided as forcing files

    subglacial_discharge = 0.0d0

    do j = 1, ny
       do i = 1, nx
          nb = basin_number(i,j)
          if (nb >= 1) then
             if (area_submerged_sum_basin(nb) > 0.0d0) then
                ! Divide basin runoff (m^3/s) by basin-wide submerged area (m^2)
                subglacial_discharge(i,j) = runoff_sum_basin(nb) / area_submerged_sum_basin(nb)  ! m/s
             endif
          endif
       enddo   ! i
    enddo   ! j

  end subroutine glissade_subglacial_discharge

!-------------------------------------------------------------------------------

  subroutine glissade_thermal_forcing_avg_3d_to_2d(&
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

    integer, dimension(0:nzocn) :: zbnd   ! depths of layer boundaries

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

    thermal_forcing_2d(i,j) = 0.0d0

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
                if (thermal_forcing(k,i,j) > -99998) then    !TODO - Rewrite
                   !WHL - Limit TF to be non-negative; is that correct?
                   thermal_forcing_layer = max(thermal_forcing(k,i,j), 0.0d0)
                   thermal_forcing_2d(i,j) = thermal_forcing_2d(i,j) + thermal_forcing_layer*dlayer*layer_frac
                endif
             enddo
          enddo
       endif
    enddo

    ! Divide by the depth range
    where (thermal_forcing_2d > 0.0d0)
       thermal_forcing_2d = thermal_forcing_2d / (ztop - zbot)
    endwhere

  end subroutine glissade_thermal_forcing_avg_3d_to_2d

!-------------------------------------------------------------------------------

end module glissade_lateral_melt

!-------------------------------------------------------------------------------
