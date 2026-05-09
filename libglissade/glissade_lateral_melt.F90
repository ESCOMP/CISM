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
  use glimmer_physcon, only: rhoi, rhoo, grav, scyr
  use glimmer_log
  use glimmer_utils, only: point_diag

  use cism_parallel, only: this_rank, main_task, nhalo, &
       parallel_halo, parallel_globalindex
!  use cism_parallel, only:, parallel_global_sum, &
!       parallel_reduce_sum, parallel_reduce_max, parallel_reduce_log_or

  implicit none

  private
  public :: glissade_lateral_melt_solve
!  public :: average_thermal_forcing

  public :: verbose_latmelt

  logical, parameter :: verbose_latmelt = .true.

contains

!-------------------------------------------------------------------------------
!TODO - Remove this subroutine, if lateral melt will be subsumed under the calving solve
  
  subroutine glissade_lateral_melt_solve(model)

    !HG: adding a fullgrid submarine melt parameterisation for Greenland marine-terminated margins.
    !  It operates similar to the calving process.
    !  The cases implemented below (exept for LATERAL_MELT_NONE) should be used with CALVING_FLOAT_ZERO.

    use glissade_masks, only : glissade_get_masks

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! local variables

    integer :: nx, ny               ! horizontal grid dimensions 
    integer :: itest, jtest, rtest  ! coordinates of diagnostic point 
    real(dp) :: dx, dy              ! cell dimensions in x and y directions (m)
    real(dp) :: dt                  ! timestep (s)
    real(dp) :: time                ! current time (yr)

    type(parallel_type) :: parallel   ! info for parallel communication

    type(glide_lateral_melt) :: lateral_melt
    
    ! basic masks
    integer, dimension(model%general%ewn, model%general%nsn)  ::  &
         ice_mask,               & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,          & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,             & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask                 ! = 1 where topg is at or above sea level, else = 0

    ! subgrid masks
    integer, dimension(model%general%ewn, model%general%nsn)  ::  &
         partial_mf_mask,      & ! = 1 for partially filled MF cells (thck < thck_effective), else = 0
         full_mask               ! = 1 for ice-filled cells that are not partial_mf cells, else = 0 

    real(dp), dimension(model%general%ewn, model%general%nsn)  :: &
         mf_length               ! length of melt front within a cell

    ! Initialize
    nx = model%general%ewn
    ny = model%general%nsn
    dx = model%numerics%dew
    dy = model%numerics%dns
    dt = model%numerics%dt
    time = model%numerics%time
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    parallel = model%parallel
    lateral_melt = model%lateral_melt

    if (model%options%which_lateral_melt == LATERAL_MELT_NONE) then
       model%lateral_melt%melt_thck = 0.0d0
       if (verbose_latmelt .and. main_task) write(iulog,*) 'No lateral melt at cliff fronts'
       return
    endif

    ! Prep for other lateral melt cases

    call parallel_halo(model%geometry%thck, parallel)

    ! Get masks.
    ! Use thickness limit of 0.0 instead of thklim so as to apply to ice from any cell
    !   not just dynamically active ice.       
    call glissade_get_masks(&
         nx,            ny,             &
         parallel,                      &
         model%geometry%thck,           &
         model%geometry%topg,           &
         model%climate%eus,             &
         0.0d0,                         &   ! thklim = 0.0  !TODO - eps11?
         ice_mask,                      &
         floating_mask = floating_mask, &
         ocean_mask = ocean_mask,       &
         land_mask = land_mask)

    call parallel_halo(ocean_mask, parallel)

    select case(model%options%which_lateral_melt)

    case(LATERAL_MELT_CONSTANT)

    case(LATERAL_MELT_ISMIP6)

    case(LATERAL_MELT_COUPLED)

    end select

  end subroutine glissade_lateral_melt_solve

!-------------------------------------------------------------------------------

  !TODO - Call this subroutine from glissade_calving
  !       Pass in cf_length and return melt_thck

  subroutine constant_lateral_melt(&
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
       melt_thck)                           ! m

    ! Apply lateral melt "horizontally" based on a prescribed constant melt rate

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

    real(dp), dimension(nx,ny), intent(in) :: & 
         topg,                   & ! bedrock elevation (m)
         mf_length                 ! length of melt front in each grid cell (m)

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck                      ! ice thickness (m)

    real(dp), intent(in) :: eus    ! eustatic sea level (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         melt_thck                 ! thickness reduction (m) due to lateral melt

    ! local variables

    integer :: i, j

    !WHL - Is m_sr needed?
    real(dp) :: &
         m_sr                      ! horizontal melting rate in m/yr calculated from Slater ISMIP6 melt approach

    !TODO - Pass this in?
    real(dp), dimension(nx,ny) :: &
         thck_effective            ! thickness (m) of submerged ice

    ! Initialize

    melt_thck = 0.0d0

    !WHL - commented out the following, since assuming the melt front is grounded for now. 
    !TODO - Modify to allow floating ice at the melt front?
    !WHL - Need to check this with Heiko
!!    ! submerged thickness: flotation thickness capped by ground below water
!!    thck_effective = min(thck*(rhoi/rhoo), max(eus-topg,0.))

    ! Compute the submerged ice thickness
    ! Set to the negative of the topography for marine-grounded ice.
    ! Set to zero for land-grounded ice.

    thck_effective = max(eus-topg, 0.0d0)

    ! Loop over locally owned cells
    ! Melt occurs only in MF cells: marine ice-filled cells with one or more ocean neighbors.
    
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (melt_front_mask(i,j) == 1) then

             ! Apply an absolute “horizontal” melt rate - i.e., the melt rate perpendicular
             ! to the face of an approximately vertical calving front
!!             m_sr = frontal_melt_rate ! m/s
!!             melt_thck(i,j) = min((m_sr*dt * thck_effective(i,j) * cf_length(i,j)) / (dx*dy), thck(i,j))

             !WHL - commented out m_sr; just using the prescribed melt rate with thck_submerged
             !TODO - Question: What happens to the ice above sea level? Should it collapse?
             melt_thck(i,j) = melt_rate_const*dt * thck_effective(i,j) * mf_length(i,j) / (dx*dy)

             !TODO - Modify to extend upstream if all the ice melts
             melt_thck(i,j) = min(melt_thck(i,j), thck(i,j))
             thck(i,j) = thck(i,j) - melt_thck(i,j)
             
             if (verbose_latmelt) then
                if (this_rank == rtest .and. i == itest .and. j == jtest) then
                   write(iulog,*) 'Constant lateral melting, rank, i, j:', this_rank, i, j
                   write(iulog,*) 'H, H_eff, topg:', thck(i,j), thck_effective(i,j), topg(i,j)-eus
                   write(iulog,*) 'rate (m/yr), mf_length, melt_thck:', melt_rate_const*scyr, mf_length(i,j), melt_thck(i,j)
                endif
                call point_diag(melt_thck, 'lateral melt_thck', itest, jtest, rtest, 7, 7)
             endif

          endif   ! melt_front_mask = 1
       enddo   ! i
    enddo   ! j

  end subroutine constant_lateral_melt

!-------------------------------------------------------------------------------

  subroutine ismip6_lateral_melt(&
       nx,                 ny,           &
       dx,                 dy,           &
       dt,                 time,         &  ! s
       itest,   jtest,     rtest,        &
       melt_factor,                      &
       melt_front_mask,                  &
       subglacial_discharge,             &  ! m/s
       tforcing_2d,                      &  ! K
       thck,                             &  ! m
       topg,                             &  ! m
       eus,                              &  ! m
       mf_length,                        &  ! m
       melt_thck)                           ! m

    ! Apply lateral melt horizontally as a function of subglacial disharge and thermal forcing.
    ! Based on the parameterization of X.

    ! input/output arguments
    
    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt,                     & ! time step (s)
         time                      ! elapsed time (s) of model run

    real(dp), intent(in) :: &
         melt_factor               ! multiplier for Rignot melt parameterisation

    integer, dimension(nx,ny), intent(in)  ::  &
!!         melt_front_mask           ! = 1 where ice is grounded below sea level or floating  !HG version
         melt_front_mask           ! = 1 where ice is grounded below sea level
                                   ! and borders at least one ocean cell, else = 0

    !WHL - How is discharge computed with units of m/s?
    real(dp), dimension(nx,ny), intent(in) :: & 
         subglacial_discharge,   & ! subglacial meltwater discharge (m/s)
         tforcing_2d,            & ! average thermal forcing over some depth range
         topg,                   & ! bedrock elevation (m)
         mf_length                 ! length of melt front in each grid cell (m)

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck                      ! ice thickness (m)

    real(dp), intent(in) :: eus    ! eustatic sea level (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         melt_thck                 ! thickness reduction (m) due to lateral melt

    ! local variables

    integer :: i, j

    real(dp) :: &
         q_sr,                   & ! runoff in Rignot calculation in m/d.  
         tf_sr,                  & ! thermal forcing in deg C
         m_sr                      ! melting rate in m/yr calculated from Slater ISMIP6 melt approach

    !TODO - Pass this in?
    real(dp), dimension(nx,ny) :: &
         thck_effective            ! effective thickness (m) of submerged ice

        ! Initialize

    melt_thck = 0.0d0

    ! Compute the submerged ice thickness
    ! Set to the negative of the topography for marine-grounded ice.
    ! Set to zero for land-grounded ice.
    !TODO - Allow a floating margin

    thck_effective = max(eus-topg, 0.0d0)

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
             q_sr = subglacial_discharge(i,j) * 86400. ! runoff_applied passed in m/s; for Rignot equation convert to [m/d]

             ! Rignot et al. 2016; formulted in m/d, converted to m/s. Mulitplier frontal_melt_factor as proposed for ISMIP7
             m_sr = melt_factor * (3.0d-4 * thck_effective(i,j) * q_sr**0.39 + 0.15) * tf_sr**1.18 * 365./scyr 

             ! calculate applied thickness change
             melt_thck(i,j) = m_sr*dt * thck_effective(i,j) * mf_length(i,j) / (dx*dy)

             ! limit by local thickness
             !TODO - Do not limit; allow melting to continue upstream
             melt_thck(i,j) = min(melt_thck(i,j), thck(i,j))

             ! Update thickness 
             !TODO - Change the thickness later, in the calving calculation
             thck(i,j) = thck(i,j) - melt_thck(i,j)

             if (verbose_latmelt) then
                if (this_rank == rtest .and. i == itest .and. j == jtest) then
                   write(iulog,*) 'ISMIP6 lateral melting, rank, i, j:', this_rank, i, j
                   write(iulog,*) 'H, H_eff, topg:', thck(i,j), thck_effective(i,j), topg(i,j)-eus
                   write(iulog,*) 'mf_length, melt_thck:', mf_length(i,j), melt_thck(i,j)
                endif
                call point_diag(melt_thck, 'lateral melt_thck', itest, jtest, rtest, 7, 7)
             endif

          endif
       enddo
    enddo

  end subroutine ismip6_lateral_melt

!-------------------------------------------------------------------------------

  subroutine average_thermal_forcing

  end subroutine average_thermal_forcing

!-------------------------------------------------------------------------------

end module glissade_lateral_melt

!-------------------------------------------------------------------------------
