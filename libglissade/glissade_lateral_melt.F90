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
  public :: glissade_lateral_melt_constant, glissade_lateral_melt_ismip6, &
       glissade_lateral_thermal_forcing_avg

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
       thck,                             &  ! m
       topg,                             &  ! m
       eus,                              &  ! m
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
!!         melt_front_mask           ! = 1 where ice is grounded below sea level or floating  !HG version
         melt_front_mask           ! = 1 where ice is grounded below sea level
                                   ! and borders at least one ocean cell, else = 0

    real(dp), intent(in) :: &
         melt_rate_const           ! prescribed constant melt rate (m/yr)

    !Note: Typically, this is the same as the calving front length
    real(dp), dimension(nx,ny), intent(in) :: & 
         topg,                   & ! bedrock elevation (m)
         mf_length                 ! length of melt front in each grid cell (m)

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck                      ! ice thickness (m); typically = thck_effective from subgrid CF scheme

    real(dp), intent(in) :: eus    ! eustatic sea level (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         latmelt_dthck             ! thickness reduction (m) due to lateral melt

    ! local variables

    integer :: i, j

!!    real(dp) :: &
!!         m_sr                      ! horizontal melting rate in m/yr calculated from Slater ISMIP6 melt approach

    real(dp), dimension(nx,ny) :: &
         thck_submerged            ! effective thickness (m) of submerged ice

    ! Initialize

    latmelt_dthck = 0.0d0

    ! Compute the submerged ice thickness
    ! Set to the negative of the topography for marine-grounded ice.
    ! Set to zero for land-grounded ice.

    thck_submerged = thck*(rhoi/rhoo)
    thck_submerged = min(thck_submerged, max(eus-topg,0.0d0))

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

             !TODO - Modify to extend upstream if all the ice melts
!             latmelt_dthck(i,j) = min(latmelt_dthck(i,j), thck(i,j))
!             thck(i,j) = thck(i,j) - latmelt_dthck(i,j)
             
             if (verbose_latmelt) then
                if (this_rank == rtest .and. i == itest .and. j == jtest) then
                   write(iulog,*) 'Constant lateral melting, rank, i, j:', this_rank, i, j
                   write(iulog,*) 'H, H_eff, topg:', thck(i,j), thck_submerged(i,j), topg(i,j)-eus
                   write(iulog,*) 'rate (m/yr), mf_length, latmelt_dthck:', &
                        melt_rate_const*scyr, mf_length(i,j), latmelt_dthck(i,j)
                endif
                call point_diag(latmelt_dthck, 'lateral melt dthck', itest, jtest, rtest, 7, 7)
             endif

          endif   ! melt_front_mask = 1
       enddo   ! i
    enddo   ! j

  end subroutine glissade_lateral_melt_constant

!-------------------------------------------------------------------------------

  subroutine glissade_lateral_melt_ismip6(&
       nx,                 ny,           &
       dx,                 dy,           &
       dt,                 time,         &  ! s
       itest,   jtest,     rtest,        &
       melt_front_mask,                  &
       melt_factor,                      &
       subglacial_discharge,             &  ! m/s
       tforcing_2d,                      &  ! K
       thck,                             &  ! m
       topg,                             &  ! m
       eus,                              &  ! m
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

    !WHL - How is discharge computed with units of m/s?
    real(dp), dimension(nx,ny), intent(in) :: & 
         subglacial_discharge,   & ! subglacial meltwater discharge (m/s)
         tforcing_2d,            & ! average thermal forcing over some depth range (K)
         topg,                   & ! bedrock elevation (m)
         mf_length                 ! length of melt front in each grid cell (m)

    real(dp), dimension(nx,ny), intent(in) :: &
         thck                      ! ice thickness (m); typically = thck_effective from subgrid CF scheme

    real(dp), intent(in) :: eus    ! eustatic sea level (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         latmelt_dthck             ! thickness reduction (m) due to lateral melt

    ! local variables

    integer :: i, j

    real(dp) :: &
         q_sr,                   & ! runoff in Rignot calculation in m/d
         tf_sr,                  & ! thermal forcing in deg C
         m_sr                      ! melting rate in m/yr calculated from Slater ISMIP6 melt approach

    !TODO - Pass this in?
    real(dp), dimension(nx,ny) :: &
         thck_submerged            ! effective thickness (m) of submerged ice

        ! Initialize

    latmelt_dthck = 0.0d0

    ! Compute the submerged ice thickness
    ! Set to the negative of the topography for marine-grounded ice.
    ! Set to zero for land-grounded ice.

    thck_submerged = thck*(rhoi/rhoo)
    thck_submerged = min(thck_submerged, max(eus-topg,0.0d0))

    ! Loop over locally owned cells
    ! Melt occurs only in MF cells: marine ice-filled cells with one or more ocean neighbors.
    
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (melt_front_mask(i,j) == 1) then

             ! Apply an absolute “horizontal” melt rate - i.e., the melt rate perpendicular
             !  to the face of an approximately vertical calving front
             ! Parameterize melt as a function of subglacial discharge, thermal forcing and effective thickness
             ! See Rignot et al. 2016 or ISMIP6 melt forcing approach.

             tf_sr = tforcing_2d(i,j) ! 2d thermal forcing [degC]
             q_sr = subglacial_discharge(i,j) * scday ! runoff_applied passed in m/s; for Rignot equation convert to [m/d]

             ! Rignot et al. 2016; formulated in m/d, converted to m/s. Mulitplier frontal_melt_factor as proposed for ISMIP7
             m_sr = melt_factor * (3.0d-4 * thck_submerged(i,j) * q_sr**0.39d0 + 0.15d0) * tf_sr**1.18d0 * 365.d0/scyr

             ! calculate applied thickness change
             latmelt_dthck(i,j) = m_sr * dt * thck_submerged(i,j) * mf_length(i,j) / (dx*dy)

             ! limit by local thickness
             !TODO - Do not limit; allow melting to continue upstream in the calving calculation
!!             latmelt_dthck(i,j) = min(latmelt_dthck(i,j), thck(i,j))

             ! Update thickness 
             !WHL - Change the thickness later, in the calving calculation
!             thck(i,j) = thck(i,j) - latmelt_dthck(i,j)

             if (verbose_latmelt) then
                if (this_rank == rtest .and. i == itest .and. j == jtest) then
                   write(iulog,*) 'ISMIP6 lateral melting, rank, i, j:', this_rank, i, j
                   write(iulog,*) 'H, H_sub, topg:', thck(i,j), thck_submerged(i,j), topg(i,j)-eus
                   write(iulog,*) 'mf_length, latmelt_dthck:', mf_length(i,j), latmelt_dthck(i,j)
                endif
                call point_diag(latmelt_dthck, 'lateral melt dthck', itest, jtest, rtest, 7, 7)
             endif

          endif
       enddo
    enddo

  end subroutine glissade_lateral_melt_ismip6

!-------------------------------------------------------------------------------

  subroutine glissade_lateral_thermal_forcing_avg(&
       nx,              ny,      &
       nzocn,                    &
       zocn,                     &
       thermal_forcing,          &
       tforcing_2d,              &
       ztop_in,         zbot_in)

    ! Average the thermal forcing over a prescribed depth range.
    ! The default range is -200 m to -500 m.

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny                    !> number of grid cells in each dimension

    integer, intent(in) :: &
         nzocn                     !> number of ocean levels

    real(dp), dimension(nzocn), intent(in) :: &
         zocn                      !> ocean levels (m) where forcing is provided, negative below sea level

    real(dp), dimension(nzocn,nx,ny), intent(in) :: &
         thermal_forcing           !> thermal forcing field at ocean levels

    real(dp), dimension(nx,ny), intent(out) :: &
         tforcing_2d               !> average thermal forcing

    real(dp), intent(in), optional :: &
         ztop_in, zbot_in          !> top and bottom of depth range (m), negative below sea level

    ! local variables

    integer :: i, j, k
    real(dp) :: ztop, zbot         ! local versions of ztop_in, zbot_in
    real(dp) :: &
         layer_frac,             & ! fraction of layer within the depth range
         dlayer,                 & ! layer thickness
         tforcing_layer            ! thermal forcing in the layer

    integer, dimension(0:nzocn) :: zbnd   ! depths of layer boundaries

    ! ISMIP values for the depth range
    real(dp), parameter :: ztop_ismip = -200.d0    ! top of depth range (m), ISMIP6 parameterization
    real(dp), parameter :: zbot_ismip = -500.d0    ! bottom of depth range (m), ISMIP6 parameterization

    if (present(ztop_in)) then
       ztop = ztop_in
    else
       ztop = ztop_ismip
    endif

    if (present(zbot_in)) then
       zbot = zbot_in
    else
       zbot = zbot_ismip
    endif

    if (ztop >= 0.0d0 .or. zbot >= 0.0d0) then
       call write_log('Error, average_thermal_forcing, zbot and ztop must be < 0', GM_FATAL)
    endif

    ! initialize
    tforcing_2d = 0.0d0

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

    tforcing_2d(i,j) = 0.0d0

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
                   tforcing_layer = max(thermal_forcing(k,i,j), 0.0d0)
                   tforcing_2d(i,j) = tforcing_2d(i,j) + tforcing_layer*dlayer*layer_frac
                endif
             enddo
          enddo
       endif
    enddo

    ! Divide by the depth range
    where (tforcing_2d > 0.0d0)
       tforcing_2d = tforcing_2d / (ztop - zbot)
    endwhere

  end subroutine glissade_lateral_thermal_forcing_avg

!-------------------------------------------------------------------------------

end module glissade_lateral_melt

!-------------------------------------------------------------------------------
