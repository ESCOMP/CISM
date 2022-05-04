!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_calving.F90 - part of the Community Ice Sheet Model (CISM)  
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
!!#ifdef HAVE_CONFIG_H
!!#include "config.inc"
!!#endif

module glissade_calving

  use glide_types
  use glimmer_global, only: dp
  use glimmer_log
  use cism_parallel, only: this_rank, main_task, nhalo, &
       parallel_halo, parallel_globalindex, parallel_reduce_sum, parallel_reduce_max

  use glimmer_paramets, only: eps08, thk0
  use glimmer_physcon, only: rhoi, rhoo, grav, scyr
  use glide_diagnostics, only: point_diag

  implicit none

  private
  public :: glissade_calving_mask_init, glissade_thck_calving_threshold_init, &
            glissade_calve_ice, glissade_remove_icebergs, glissade_remove_isthmuses, &
            glissade_cull_calving_front, glissade_limit_cliffs
  public :: verbose_calving

!!  logical, parameter :: verbose_calving = .false.
  logical, parameter :: verbose_calving = .true.

contains

!-------------------------------------------------------------------------------

  subroutine glissade_calving_mask_init(dx,                dy,               &
                                        parallel,                            &
                                        thck,              topg,             &
                                        eus,               thklim,           &
                                        usfc_obs,          vsfc_obs,         &
                                        calving_front_x,   calving_front_y,  &
                                        calving_mask)

    ! Compute an integer calving mask if needed for the CALVING_GRID_MASK option

    use glissade_masks, only: glissade_get_masks

    ! Input/output arguments

    real(dp), intent(in) :: dx, dy                 !> cell dimensions in x and y directions (m)
    type(parallel_type), intent(in) :: parallel    !> info for parallel communication
    real(dp), dimension(:,:), intent(in) :: thck   !> ice thickness (m)
    real(dp), dimension(:,:), intent(in) :: topg   !> present bedrock topography (m)
    real(dp), intent(in) :: eus                    !> eustatic sea level (m)
    real(dp), intent(in) :: thklim                 !> minimum thickness for dynamically active grounded ice (m)
    real(dp), dimension(:,:), intent(in) :: &
         usfc_obs, vsfc_obs                        !> observed surface velocity components (m/yr)
    real(dp), intent(in) :: calving_front_x        !> for CALVING_GRID_MASK option, calve ice wherever abs(x) > calving_front_x (m)
    real(dp), intent(in) :: calving_front_y        !> for CALVING_GRID_MASK option, calve ice wherever abs(y) > calving_front_y (m)

    integer, dimension(:,:), intent(inout) :: calving_mask   !> output mask: calve floating ice wherever the mask = 1

    ! Local variables

    real(dp) :: xcell, ycell     ! global cell center coordinates (m)
    integer :: nx, ny            ! horizontal grid dimensions
    integer :: i, j              ! local cell indices
    integer :: iglobal, jglobal  ! global cell indices

    integer, dimension(:,:), allocatable :: &
         ice_mask,             & ! = 1 where ice is present
         ocean_mask              ! = 1 for ice-free ocean

    real(dp) :: mask_maxval      ! maxval of calving_mask

    nx = size(calving_mask,1)
    ny = size(calving_mask,2)

    mask_maxval = maxval(calving_mask)
    mask_maxval = parallel_reduce_max(mask_maxval)

    ! Compute the calving mask, if not read in at initialization
 
    if (mask_maxval > 0) then

       ! calving_mask was read from the input file; do not need to compute a mask here

       if (main_task) write(6,*) 'Calving_mask was read from the input file'

    elseif (calving_front_x > 0.0d0 .or. calving_front_y > 0.0d0) then

       if (main_task) write(6,*) 'Computing calving_mask based on calving_front_x/y'

       ! initialize
       calving_mask(:,:) = 0   ! no calving by default

       if (calving_front_x > 0.0d0) then

          ! set calving_mask = 1 where abs(x) > calving_front_x

          do j = 1, ny
             do i = 1, nx

                ! find global i and j indices
                call parallel_globalindex(i, j, iglobal, jglobal, parallel)

                ! find cell center x coordinate
                xcell = (dble(iglobal) - 0.5d0) * dx

                ! set calving mask = 1 based on cell coordinates relative to the calving front
                ! Note: Using absolute value to support symmetry with respect to x = 0
                if (abs(xcell) > calving_front_x) then
                   calving_mask(i,j) = 1
                endif

             enddo   ! i
          enddo   ! j

       endif   ! calving_front_x > 0

       if (calving_front_y > 0.0d0) then

          ! set calving_mask = 1 where abs(y) > calving_front_y

          do j = 1, ny
             do i = 1, nx

                ! find global i and j indices
                call parallel_globalindex(i, j, iglobal, jglobal, parallel)

                ! find cell center y coordinate
                ycell = (dble(jglobal) - 0.5d0) * dy

                ! set calving mask = 1 based on cell coordinates relative to the calving front
                if (abs(ycell) > calving_front_y) then
                   calving_mask(i,j) = 1
                endif

             enddo   ! i
          enddo   ! j

       endif   ! calving_front_y > 0

    else  ! compute the calving mask based on the initial ice extent
 
       if (main_task) then
          write(6,*) 'Computing calving_mask based on initial ice extent'
       endif

       ! initialize
       calving_mask(:,:) = 0  ! no calving by default

       ! Get an ocean mask
       allocate(ice_mask(nx,ny))
       allocate(ocean_mask(nx,ny))

       !TODO: Modify glissade_get_masks so that 'parallel' is not needed
       call glissade_get_masks(nx,            ny,             &
                               parallel,                      &
                               thck,          topg,           &
                               eus,           thklim,         &
                               ice_mask,                      &
                               ocean_mask = ocean_mask)

       ! Set the calving mask to include all ice-free ocean cells.
       ! Make an exception for cells where usfc_obs or vsfc_obs > 0.
       ! This would include cells with observed nonzero velocity (and hence ice present)
       !  which are ice-free ocean in the input thickness dataset (e.g., Bedmachine).
       !  As of Dec. 2021, this is the case for parts of the Thwaites shelf region.
       !  We want to allow the shelf to expand into regions where ice was present
       !   and flowing recently, even if no longer present.
       ! Any ice entering these cells during the run will calve.

       do j = 2, ny-1
          do i = 2, nx-1
             if (ocean_mask(i,j) == 1) then
                if (usfc_obs(i-1,j)   == 0.0d0 .and. usfc_obs(i,j)   == 0.0d0 .and. &
                    usfc_obs(i-1,j-1) == 0.0d0 .and. usfc_obs(i,j-1) == 0.0d0 .and. &
                    vsfc_obs(i-1,j)   == 0.0d0 .and. vsfc_obs(i,j)   == 0.0d0 .and. &
                    vsfc_obs(i-1,j-1) == 0.0d0 .and. vsfc_obs(i,j-1) == 0.0d0) then
                   calving_mask(i,j) = 1   ! calve ice in this cell
                else
                   calving_mask(i,j) = 0
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   write(6,*) 'ocean cell with uobs, vobs > 0: iglobal, jglobal, thck, uobs, vobs', &
                        iglobal, jglobal, thck(i,j), usfc_obs(i,j), vsfc_obs(i,j)
                endif
             else
                calving_mask(i,j) = 0
             endif
          enddo
       enddo

       call parallel_halo(calving_mask, parallel)

       deallocate(ice_mask)
       deallocate(ocean_mask)

    endif  ! mask_maxval > 0

    ! halo update moved to higher level
    call parallel_halo(calving_mask, parallel)

  end subroutine glissade_calving_mask_init

!-------------------------------------------------------------------------------

  subroutine glissade_thck_calving_threshold_init(&
       nx,      ny,               &
       parallel,                  &
       itest,   jtest,   rtest,   &
       which_ho_calving_front,    &
       thck,    topg,             &
       eus,     thklim,           &
       marine_connection_mask,    &
       calving_minthck,           &
       thck_calving_threshold)

    ! Given the input ice thickness, identify calving-front cells and compute the CF thickness.
    ! This subroutine is called at initialization when whichcalving = CALVING_THCK_THRESHOLD.
    !
    ! There are two options:
    ! (1) Set thck_calving_threshold to a uniform value of calving_minthck.
    ! (2) Set thck_calving_threshold to thck_calving_front, as computed from the input ice thickness.
    ! Note: This subroutine uses SI units.

    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask
    use glissade_grid_operators, only: glissade_scalar_extrapolate

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

    integer, intent(in) :: nx, ny                              !> horizontal grid dimensions
    type(parallel_type), intent(in) :: parallel                !> info for parallel communication
    integer, intent(in) :: itest, jtest, rtest                 !> coordinates of diagnostic point
    integer, intent(in) :: which_ho_calving_front              !> = 1 for subgrid calving-front scheme, else = 0

    real(dp), dimension(nx,ny), intent(in) :: thck             !> ice thickness (m)
    real(dp), dimension(nx,ny), intent(in) :: topg             !> present bedrock topography (m)
    real(dp), intent(in)                   :: eus              !> eustatic sea level (m)
    real(dp), intent(in)                   :: thklim           !> minimum thickness for dynamically active grounded ice (m)
    integer,  dimension(nx,ny), intent(in) :: marine_connection_mask  !> = 1 for cells with a marine connection to the ocean

    ! Note: calving_minthck = 0 by default.
    !       If calving_minthck > 0 is specified in the config file, then we set thck_calving_threshold = calving_minthck.
    !       Otherwise, we set thck_calving_threshold to the input ice thickness (thck)

    real(dp), intent(in)                   :: calving_minthck  !> uniform thickness threshold for calving ice (m)
    real(dp), dimension(nx,ny), intent(out) ::  &
         thck_calving_threshold                 !> ice in the calving domain will calve when thinner than this value (m)

    ! Local variables

    integer :: i, j, n

    integer, dimension(nx,ny) :: &
         ice_mask,                  & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,             & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,                & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,                 & ! = 1 where topg is at or above sea level, else = 0
         calving_front_mask           ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny) :: &
         thck_calving_front,             & ! effective ice thickness at the calving front
         thck_calving_threshold_smoothed   ! smoothed version of thck_calving_threshold

    !TODO - Make the calving thresholds config parameters?
    real(dp), parameter :: &
         calving_threshold_min = 100.d0,  &! min allowed value (m) of thck_calving_threshold
         calving_threshold_max = 500.d0    ! max allowed value (m) of thck_calving_threshold

    if (calving_minthck > eps08) then

       ! Set thck_calving_threshold to a uniform value
       thck_calving_threshold(:,:) = calving_minthck

       if (verbose_calving .and. main_task) then
          write(6,*) 'Set thck_calving_threshold to calving_minthck (m):', calving_minthck
       endif

    else

       ! Set thck_calving_threshold based on the input ice thickness and calving-front thickness

       ! Get masks

       call glissade_get_masks(nx,            ny,             &
                               parallel,                      &
                               thck,          topg,           &
                               eus,           thklim,         &
                               ice_mask,                      &
                               floating_mask = floating_mask, &
                               ocean_mask = ocean_mask,       &
                               land_mask = land_mask)

       ! Note: If using a subgrid CF scheme, thck_calving_front is the effective ice thickness at the calving front,
       !        equal to the mean thickness of marine interior neighbors.
       !       Without a subgrid CF scheme, thck_calving_front is set to thck in calving_front cells.

       call glissade_calving_front_mask(nx,            ny,                 &
                                        which_ho_calving_front,            &
                                        parallel,                          &
                                        thck,          topg,               &
                                        eus,                               &
                                        ice_mask,      floating_mask,      &
                                        ocean_mask,    land_mask,          &
                                        calving_front_mask,                &
                                        thck_calving_front)

       if (verbose_calving) then
          call point_diag(thck, 'Set thck_calving_threshold, thck, (m)', itest, jtest, rtest, 7, 7)
          call point_diag(floating_mask, 'floating_mask', itest, jtest, rtest, 7, 7)
          call point_diag(calving_front_mask, 'calving_front_mask', itest, jtest, rtest, 7, 7)
          call point_diag(thck_calving_front, 'thck_calving_front (m)', itest, jtest, rtest, 7, 7)
       endif

       if (which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then

          ! Extrapolate and smooth thck_calving_front, starting with CF cells where it was just computed

          if (verbose_calving .and. main_task) then
             write(6,*) ' '
             write(6,*) 'Set thck_calving_threshold by extrapolating thck_calving_front to marine-connected cells'
          endif

          ! Limit thck_calving_front to lie within a prescribed range.
          where (calving_front_mask == 1)
             thck_calving_front = max(thck_calving_front, calving_threshold_min)
             thck_calving_front = min(thck_calving_front, calving_threshold_max)
          endwhere

          if (verbose_calving) then
             call point_diag(thck_calving_front, 'After corrections, thck_calving_front (m)', &
                  itest, jtest, rtest, 7, 7)
          endif

          ! Extrapolate the CF thickness from cells with calving_front_mask = 1
          !  to all cells with marine_connection_mask = 1.
          ! Apply a Laplacian smoother during the extrapolation.

          call glissade_scalar_extrapolate(nx,    ny,                 &
                                           parallel,                  &
                                           calving_front_mask,        &
                                           thck_calving_front,        &
                                           marine_connection_mask,    &
                                           thck_calving_threshold,    &
                                           npoints_stencil = 9,       &
                                           apply_smoother = .true.,   &
                                           itest = itest, jtest = jtest, rtest = rtest)

          call parallel_halo(thck_calving_threshold, parallel)

          if (verbose_calving) then
             if (this_rank == rtest) then
                write(6,*) 'Extrapolated thck_calving_front to interior marine-based cells'
             endif
             call point_diag(thck_calving_threshold, 'thck_calving_threshold (m)', itest, jtest, rtest, 7, 7)
          endif

       elseif (which_ho_calving_front == HO_CALVING_FRONT_NO_SUBGRID) then

          ! Set thck_calving_threshold = thck in marine-connected cells.
          ! As a result, any floating ice at the calving front, once thinned by mass-balance
          !  or ice-dynamics changes, will experience increased calving.

          where (marine_connection_mask == 1)
             thck_calving_threshold = thck
          elsewhere
             thck_calving_threshold = 0.0d0
          endwhere

          if (verbose_calving) then
             if (this_rank == rtest) then
                write(6,*) 'Set thck_calving_threshold = thck in marine-connected cells'
             endif
             call point_diag(thck_calving_threshold, 'thck_calving_threshold (m)', itest, jtest, rtest, 7, 7)
          endif

       endif   ! which_ho_calving_front

    endif  ! calving_minthck > eps08

  end subroutine glissade_thck_calving_threshold_init

!-------------------------------------------------------------------------------

  subroutine glissade_calve_ice(nx,               ny,    &
                                which_calving,           &
                                calving_domain,          &
                                which_ho_calving_front,  &
                                parallel,                &
                                calving,                 &  ! calving derived type
                                itest,   jtest,   rtest, &
                                dt,                      &  ! s
                                dx,               dy,    &  ! m
                                sigma,                   &
                                thklim,                  &  ! m
                                velnorm_mean,            &  ! m/s
                                thck,             relx,  &  ! m
                                topg,             eus)      ! m

    ! Calve ice according to one of several methods.
    ! Note: This subroutine uses SI units.

    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask

    implicit none

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

    integer, intent(in) :: nx, ny                  !> horizontal grid dimensions

    integer, intent(in) :: which_calving           !> option for calving law
    integer, intent(in) :: calving_domain          !> option for where calving can occur
                                                   !> = 0 if calving occurs at the ocean edge only
                                                   !> = 1 if calving occurs everywhere the calving criterion is met
                                                   !> = 2 if calving occurs where criterion is met and there is a connected path
                                                   !>     to the ocean through other cells where the criterion is met
    integer, intent(in) :: which_ho_calving_front  !> = 1 for subgrid calving-front scheme, else = 0

    type(parallel_type), intent(in) :: parallel    !> info for parallel communication
    type(glide_calving), intent(inout) :: calving  !> calving object

!    Note: The calving object includes the following fields and parameters used in this subroutine:
!    real(dp), intent(in)                     :: marine_limit        !> lower limit on topography elevation at marine edge before ice calves
                                                                     !> Note: marine_limit (shared by Glide) has scaled model units
!    real(dp), intent(in)                     :: calving_fraction    !> fraction of ice lost at marine edge when calving; 
                                                                     !> used with CALVING_FLOAT_FRACTION
!    real(dp), intent(in)                     :: timescale           !> timescale (s) for calving; calving_thck = thck * max(dt/timescale, 1)
                                                                     !> if timescale = 0, then calving_thck = thck
!    real(dp), intent(in)                     :: thck_calving_threshold  !> calve ice in the calving domain if thck < thck_calving_threshold (m);
                                                                         !> used with CALVING_THCK_THRESHOLD, EIGENCALVING, CALVING_DAMAGE
!    real(dp), intent(in)                     :: eigencalving_constant   !> eigencalving constant; m/s (lateral calving rate) per Pa (tensile stress)
!    real(dp), intent(in)                     :: eigen2_weight       !> weight given to tau_eigen2 relative to tau_eigen1 in tau_eff (unitless)
!    real(dp), dimension(:,:), intent(in)     :: tau_eigen1          !> first eigenvalue of 2D horizontal stress tensor (Pa)
!    real(dp), dimension(:,:), intent(in)     :: tau_eigen2          !> second eigenvalue of 2D horizontal stress tensor (Pa)
!    real(dp), dimension(:,:), intent(inout)  :: tau_eff             !> effective stress (Pa) for calving; derived from tau_eigen1/2
!    real(dp), dimension(:,:,:), intent(inout):: damage              !> 3D scalar damage parameter
!    real(dp), intent(in)                     :: damage_threshold    !> threshold value where ice is sufficiently damaged to calve
!    real(dp), intent(in)                     :: damage_constant     !> rate of change of damage (1/s) per unit stress (Pa)
!    real(dp), dimension(:,:), intent(inout)  :: lateral_rate        !> lateral calving rate (m/s) at calving front,used with EIGENCALVING
!    integer,  dimension(:,:), intent(in)     :: calving_mask        !> integer mask: calve ice where calving_mask = 1
!    real(dp), dimension(:,:), intent(out)    :: calving_thck        !> thickness lost due to calving in each grid cell (m)

    integer, intent(in) :: itest, jtest, rtest                   !> coordinates of diagnostic point
    real(dp), intent(in)                      :: dt                !> model timestep (s)
    real(dp), intent(in)                      :: dx, dy            !> grid cell size in x and y directions (m)
    real(dp), dimension(:), intent(in)        :: sigma             !> vertical sigma coordinate
    real(dp), dimension(nx-1,ny-1), intent(in):: velnorm_mean      !> mean ice speed at vertices (m/s)
    real(dp), dimension(nx,ny), intent(inout) :: thck              !> ice thickness (m)
    real(dp), dimension(nx,ny), intent(in)    :: relx              !> relaxed bedrock topography (m)
    real(dp), dimension(nx,ny), intent(in)    :: topg              !> present bedrock topography (m)
    real(dp), intent(in)                      :: thklim            !> minimum thickness for dynamically active grounded ice (m)
    real(dp), intent(in)                      :: eus               !> eustatic sea level (m)

    ! local variables

    integer :: nz          ! number of vertical levels
                           ! Note: number of ice layers = nz-1
    integer :: i, j, k, n
    integer :: ii, jj

    real(dp), dimension(nx,ny) ::  &
         thck_calving_front,     & ! effective ice thickness at the calving front
         tau1, tau2,             & ! tau_eigen1 and tau_eigen2 (Pa), modified for calving
         eps1, eps2                ! eps_eigen1 and eps_eigen2 (1/s), modified for calving

    ! basic masks
    integer, dimension(nx,ny)  ::  &
         ice_mask,               & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,          & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,             & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,              & ! = 1 where topg is at or above sea level, else = 0
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    ! Note: Calving occurs in a cell if and only if (1) the calving law permits calving, 
    !       and (2) the cell is in the calving domain, as specified by the calving_domain option.
    !       The calving domain by default is limited to the ocean edge (CALVING_DOMAIN_OCEAN_EDGE), 
    !       but can be extended to include all ice-covered cells (CALVING_DOMAIN_EVERYWHERE).

    !TODO - Make these integer masks like the ones above?
    logical, dimension(nx,ny) ::  &
         calving_law_mask,    & ! = T where the calving law permits calving, else = F
         calving_domain_mask    ! = T in the domain where calving is allowed to occur (e.g., at ocean edge), else = F

    real(dp) :: &
         float_fraction_calve, & ! = calving_fraction for which_calving = CALVING_FLOAT_FRACTION
                                 ! = 1.0 for which_calving = CALVING_FLOAT_ZERO
         thinning_rate,        & ! vertical thinning rate (m/s)
         dthck                   ! thickness change (m)

    real(dp), dimension(nx,ny) :: &
         thck_effective          ! effective thickness (m) for calving, weighted toward upstream thickness

    integer, dimension(nx,ny) :: &
         partial_cf_mask,      & ! = 1 for partially filled CF cells (thck < thck_effective), else = 0
         full_mask               ! = 1 for ice-filled cells that are not partial_cf cells, else = 0

    character(len=100) :: message
   
    ! initialize

    nz = size(sigma)

    if (which_calving == CALVING_NONE) then   ! do nothing
       if (verbose_calving .and. main_task) write(6,*) 'No calving'
       return
    endif

    if (verbose_calving .and. main_task) then
       write(6,*) ' '
       write(6,*) 'In glissade_calve_ice, which_calving =', which_calving
       write(6,*) 'calving_domain =', calving_domain
    endif

    ! Set the thickness fraction to be removed in each calving cell
    ! Note: The CALVING_FLOAT_FRACTION option has been superseded by the calving%timescale variable,
    !       but is included here for consistency with Glide.
    ! TODO: Remove CALVING_FLOAT_FRACTION option?

    if (which_calving == CALVING_FLOAT_FRACTION) then

       !WHL - Changed definition of calving fraction; now it is the fraction lost
       !      rather than the fraction remaining
       float_fraction_calve = calving%calving_fraction
       
    else  ! other calving options

       if (calving%timescale == 0.0d0) then  ! calve the entire column for eligible columns (this is the default)
          float_fraction_calve = 1.0d0
       else  ! calve a fraction of the column based on the calving time scale
          float_fraction_calve = min(dt/calving%timescale, 1.0d0)
       endif
       
    endif
       
    ! Do the calving based on the value of which_calving

    ! Notes on eigencalving and damage-based calving:
    ! For both damage-based calving and eigencalving, the calving rate depends on eigenvalues of the 2D stress tensor.
    ! The main difference is that for eigencalving, the calving rate is based on current stresses
    !  at the calving front, whereas for damage-based calving, the calving rate is based on cumulative damage,
    !  which accumulates in floating cells due to stresses and then is advected downstream.

    if (which_calving == CALVING_DAMAGE) then

       call damage_based_calving(&
            nx,       ny,       nz,                    &
            dx,                 dy,                    &
            sigma,              dt,                    &
            itest,    jtest,    rtest,                 &
            parallel,                                  &
            which_ho_calving_front,                    &
            velnorm_mean,                              &  ! m/s
            thck,               topg,                  &  ! m
            eus,                thklim,                &  ! m
            calving%tau_eigen1, calving%tau_eigen2,    &  ! Pa
            calving%eps_eigen1, calving%eps_eigen2,    &  ! 1/s
            calving%damage_constant1,                  &  ! 1/s
            calving%damage_constant2,                  &  ! 1/s
            calving%damage_threshold,                  &  !
            calving%minthck,                           &  ! m
            calving%damage,                            &  !
            calving%lateral_rate,                      &  ! m/s
            calving%calving_thck)                         ! m

    elseif (which_calving == EIGENCALVING) then

       call eigen_calving(&
            nx,                 ny,                    &
            dx,                 dy,                    &
            dt,                                        &
            itest,   jtest,     rtest,                 &
            parallel,                                  &
            which_ho_calving_front,                    &
            thck,               topg,                  &  ! m
            eus,                thklim,                &  ! m
            calving%tau_eigen1, calving%tau_eigen2,    &  ! Pa
            calving%eps_eigen1, calving%eps_eigen2,    &  ! 1/s
            calving%eigenconstant1,                    &
            calving%eigenconstant2,                    &
            calving%minthck,                           &  ! m
            calving%lateral_rate,                      &  ! m/s
            calving%calving_thck)                         ! m

    elseif (which_calving == CALVING_THCK_THRESHOLD) then

       !TODO - Move thickness-based calving to subroutine

       ! Note: If running without the subgrid CF scheme, then glissade_calving_front_mask will return CF_mask = 1
       !       for floating cells adjacent to ice-free ocean, with thck_calving_front = thck.

       ! Get masks

       call glissade_get_masks(nx,            ny,             &
                               parallel,                      &
                               thck,          topg,           &
                               eus,           thklim,         &
                               ice_mask,                      &
                               floating_mask = floating_mask, &
                               ocean_mask = ocean_mask,       &
                               land_mask = land_mask)

       call glissade_calving_front_mask(&
                               nx,            ny,                 &
                               which_ho_calving_front,            &
                               parallel,                          &
                               thck,          topg,               &
                               eus,                               &
                               ice_mask,      floating_mask,      &
                               ocean_mask,    land_mask,          &
                               calving_front_mask,                &
                               thck_calving_front)

       if (verbose_calving) then
          if (this_rank == rtest) then
             write(6,*) ' '
             write(6,*) 'Thickness-based calving:'
          endif
          call point_diag(calving_front_mask, 'calving_front_mask', itest, jtest, rtest, 7, 7)
          call point_diag(thck_calving_front, 'thck_calving_front (m)', itest, jtest, rtest, 7, 7)
          call point_diag(calving%thck_calving_threshold, 'thck_calving_threshold (m)', &
               itest, jtest, rtest, 7, 7)
          call point_diag(thck, 'thck (m)', itest, jtest, rtest, 7, 7)
       endif

       ! Apply thinning in calving-front cells whose effective thickness (H_e = thck_calving_front) is less than
       !  a prescribed minimum value (Hc_min = thck_calving_threshold).
       !
       ! The effective thinning rate is given by
       !
       !    dH_e/dt = -(Hc_min - H_e) / tau_c  where Hc_min > H_e
       !    dH_e/dt = 0 elsewhere
       !
       ! where tau_c = calving%timescale.
       !
       ! The thinning rate applied to the mean cell thickness (thck) is given by
       !
       !    dH/dt = min(H/H_e, 1) * dH_e/dt
       !
       ! Thus, any ice with H_e < Hc_min is removed on a time scale given by tau_c.

       if (calving%timescale <= 0.0d0) then
          write(message,*) 'Must set calving timescale to a positive nonzero value for this calving option'
          call write_log(message, GM_FATAL)
       endif

       do j = 2, ny-1
          do i = 2, nx-1
             if (calving_front_mask(i,j) == 1 .and. &
                  thck_calving_front(i,j) > 0.0d0 .and. thck_calving_front(i,j) <= calving%thck_calving_threshold(i,j)) then
                
!!                if (verbose_calving .and. thck(i,j) > 0.0d0) &
!!                     write(6,*) 'Calve thin floating ice: task, i, j, thck =', this_rank, i, j, thck(i,j)

                ! calving%timescale has units of s
                thinning_rate = (calving%thck_calving_threshold(i,j) - thck_calving_front(i,j)) / calving%timescale
                dthck = thinning_rate * dt

                !WHL - debug
                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                   write(6,*) ' '
                   write(6,*) 'Thinning: r, i, j =', rtest, itest, jtest
                   write(6,*) 'thck:', thck(i,j)
                   write(6,*) 'thck_calving_front (m) =', thck_calving_front(i,j)
                   write(6,*) 'thck_calving_threshold (m) =', calving%thck_calving_threshold(i,j)
                   write(6,*) 'thinning rate (m/yr) =', thinning_rate * scyr
                   write(6,*) 'dthck (m) =', dthck
                endif

                if (dthck > thck(i,j)) then
                   calving%calving_thck(i,j) = calving%calving_thck(i,j) + thck(i,j)
                   thck(i,j) = 0.0d0
                else
                   thck(i,j) = thck(i,j) - dthck
                   calving%calving_thck(i,j) = calving%calving_thck(i,j) + dthck
                endif

             endif   ! thck_calving_front < thck_calving_threshold in calving_front cell
          enddo   ! i
       enddo   ! j

    else   ! other calving options
           !TODO - Put these in a separate subroutine

       ! Get masks.
       ! Use thickness limit of 0.0 instead of thklim so as to remove ice from any cell
       !  that meets the calving criteria, not just dynamically active ice.

       call glissade_get_masks(nx,            ny,             &
                               parallel,                      &
                               thck,          topg,           &
                               eus,           0.0d0,          &   ! thklim = 0.0
                               ice_mask,                      &
                               floating_mask = floating_mask, &
                               ocean_mask = ocean_mask)

       ! set the calving-law mask
       ! Note: Cells that meet the calving-law criteria will be calved provided they also lie in the calving domain,
       !       as determined below.

       select case (which_calving)

       case(CALVING_FLOAT_ZERO, CALVING_FLOAT_FRACTION)     ! calve ice that is floating

          do j = 1, ny
             do i = 1, nx
                if (floating_mask(i,j) == 1) then
                   calving_law_mask(i,j) = .true.
                else
                   calving_law_mask(i,j) = .false.
                endif
             enddo
          enddo

          !NOTE: The Glide version of CALVING_FLOAT_ZERO calves all floating ice.
          !      Glissade calves floating ice only in the calving domain, which is CALVING_DOMAIN_OCEAN_EDGE by default.
          !      Must set calving_domain = CALVING_DOMAIN_EVERYWHERE to match the Glide behavior.
          !TODO: Change the default to calving_domain_everywhere?

       case(CALVING_RELX_THRESHOLD)   ! set thickness to zero if relaxed bedrock is below a given level

          !WHL - The Glide version of CALVING_RELX_THRESHOLD calves ice wherever the relaxed bedrock criterion is met.
          !      Must set calving_domain = CALVING_DOMAIN_EVERYWHERE to match the Glide behavior.
          ! Note: calving%marine_limit (a holdover from Glide) has scaled model units
          where (relx <= calving%marine_limit*thk0 + eus)   ! convert marine_limit from scaled units to m
             calving_law_mask = .true.
          elsewhere
             calving_law_mask = .false.
          endwhere

       case(CALVING_TOPG_THRESHOLD)   ! set thickness to zero if present bedrock is below a given level

          where (topg < calving%marine_limit*thk0 + eus)    ! convert marine_limit from scaled units to m
             calving_law_mask = .true.
          elsewhere
             calving_law_mask = .false.
          endwhere

       case(CALVING_HUYBRECHTS)    ! Huybrechts grounding line scheme for Greenland initialization

          if (eus > -80.d0) then
             where (relx <= 2.d0*eus)
                calving_law_mask = .true.
             elsewhere
                calving_law_mask = .false.
             end where
          elseif (eus <= -80.d0) then
             where (relx <= (2.d0*eus - 0.25d0*(eus + 80.d0)**2.d0))
                calving_law_mask = .true.
             elsewhere
                calving_law_mask = .false.
             end where
          end if

       end select

       ! halo update (may not be necessary if thck, damage, etc. are correct in halos, but including to be safe)
       call parallel_halo(calving_law_mask, parallel)

       ! set the calving domain mask

       if (calving_domain == CALVING_DOMAIN_OCEAN_EDGE) then  ! calving domain includes floating cells at margin only
                                                              !WHL - Could modify to include grounded marine cells at margin
          do j = 2, ny-1
             do i = 2, nx-1

                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                   write(6,*) 'task, i, j, ice_mask, floating_mask:',  &
                        this_rank, i, j, ice_mask(i,j), floating_mask(i,j)
                endif

                if ( floating_mask(i,j) == 1 .and.   &
                     (ocean_mask(i-1,j)==1 .or. ocean_mask(i+1,j)==1 .or. ocean_mask(i,j-1)==1 .or. ocean_mask(i,j+1)==1) ) then
                   calving_domain_mask(i,j) = .true.
                else
                   calving_domain_mask(i,j) = .false.
                endif
             enddo
          enddo

          ! halo update (since the loop above misses some halo cells)
          call parallel_halo(calving_domain_mask, parallel)

          if (verbose_calving) then
             call point_diag(calving_domain_mask, 'calving_domain_mask', itest, jtest, rtest, 7, 7)
          endif

       elseif (calving_domain == CALVING_DOMAIN_EVERYWHERE) then  ! calving domain includes all cells

          calving_domain_mask(:,:) = .true.

       endif   ! calving_domain

       ! Calve ice where calving_law_mask = T and calving_domain_mask = T
       do j = 1, ny
          do i = 1, nx
             if (calving_law_mask(i,j) .and. calving_domain_mask(i,j)) then

                if (verbose_calving .and. this_rank==rtest .and. thck(i,j) > 0.0d0) then
!!                   write(6,*) 'Calve ice: task, i, j, calving_thck =', this_rank, i, j, float_fraction_calve * thck(i,j)
                endif

                calving%calving_thck(i,j) = calving%calving_thck(i,j) + float_fraction_calve * thck(i,j)
                thck(i,j) = thck(i,j) - float_fraction_calve * thck(i,j)
            endif
          enddo
       enddo

    endif   ! which_calving

    if (verbose_calving) then
       call point_diag(calving%calving_thck, 'calving_thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(thck, 'After calving, new thck (m)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine glissade_calve_ice

!---------------------------------------------------------------------------

  subroutine eigen_calving(&
       nx,            ny,                 &
       dx,            dy,                 &
       dt,                                &
       itest, jtest,  rtest,              &
       parallel,                          &
       which_ho_calving_front,            &
       thck,          topg,               &
       eus,           thklim,             &
       tau_eigen1,    tau_eigen2,         &
       eps_eigen1,    eps_eigen2,         &
       eigenconstant1,                    &
       eigenconstant2,                    &
       calving_minthck,                   &
       lateral_rate,                      &
       calving_thck)

    ! Compute calving based on the eigenvalues of the 2D horizontal stress tensor near the calving front.
    ! This is similar to damage-based calving, except that instead of allowing damage to accumulate,
    !  we prescribe a lateral calving rate proportional to one or more eigenvalues.
    ! The lateral rate is converted to a thinning rate based on the effective calving thickness.
    ! TODO - Computes calving rates based on stress eigenvalues only.  Add strain-rate option.

    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt                        ! time step (s)

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    integer, intent(in) :: &
         which_ho_calving_front    ! = 1 for subgrid calving-front scheme, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         topg,                   & ! bed topography (m)
         tau_eigen1, tau_eigen2, & ! eigenvalues of the horizontal stress tensor (Pa)
         eps_eigen1, eps_eigen2    ! eigenvalues of the horizontal strain rate tensor (1/s)

    real(dp), intent(in) ::  &
         eus,                    & ! eustatic sea level (m)
         thklim                    ! minimum thickness for dynamically active grounded ice (m)

    real(dp), intent(in) :: &
         eigenconstant1,         & ! lateral calving rate proportional to tau1 (m/s)
         eigenconstant2,         & ! lateral calving rate proportional to tau1 (m/s)
         calving_minthck           ! min effective thickness (m) for CF cells

    real(dp), dimension(nx,ny), intent(out) :: &
         lateral_rate              ! lateral rate of calving (m/s)

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m)
         calving_thck              ! thickness lost due to calving (m)

    ! local variables

    integer :: i, j, k, ii, jj

    integer, dimension(nx,ny)  ::  &
         ice_mask,               & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,          & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,             & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,              & ! = 1 where topg is at or above sea level, else = 0
         calving_front_mask,     & ! = 1 where ice is floating and borders at least one ocean cell, else = 0
         calving_front_mask_old

    integer, dimension(nx,ny) :: &
         partial_cf_mask,        & ! = 1 for partially filled CF cells (thck < thck_effective), else = 0
         full_mask                 ! = 1 for ice-filled cells that are not partial_cf cells, else = 0

    real(dp), dimension(nx,ny) :: &
         thck_calving_front,     & ! effective ice thickness at the calving front (m)
         thck_effective,         & ! effective thickness for calving (m)
         thck_old,               & ! old value of thck (m)
         dt_calving,             & ! time remaining for calving (s)
         dt_calving_old            ! old value of dt_calving

    real(dp) :: &
         ocean_neighbor_count,   & ! number of edge neighbors that are ice-free ocean
         lateral_rate_factor,    & ! multiplying factor to account for a CF longer than one edge
         thinning_rate,          & ! vertical thinning rate (m/s)
         dthck                     ! thickness change (m)

    real(dp), parameter :: &
         stress_scale_thck = 100.d0     ! ice thickness (m) in stress_scale = rhoi*g*H


    logical :: iterate_calving          ! Iterate the calving until dt_calving = 0 for all cells
    integer :: count

    ! Not sure these parallel updates are needed here
    call parallel_halo(tau_eigen1, parallel)
    call parallel_halo(tau_eigen2, parallel)
    call parallel_halo(eps_eigen1, parallel)
    call parallel_halo(eps_eigen2, parallel)

    ! Compute masks needed for calving
    ! Need a calving_front_mask; calving/thinning is applied only to cells at the calving front.

    call glissade_get_masks(&
         nx,            ny,             &
         parallel,                      &
         thck,          topg,           &
         eus,           thklim,         &
         ice_mask,                      &
         floating_mask = floating_mask, &
         ocean_mask = ocean_mask,       &
         land_mask = land_mask)

    call glissade_calving_front_mask(&
         nx,            ny,                 &
         which_ho_calving_front,            &
         parallel,                          &
         thck,          topg,               &
         eus,                               &
         ice_mask,      floating_mask,      &
         ocean_mask,    land_mask,          &
         calving_front_mask,                &
         thck_calving_front,                &
         calving_minthck = calving_minthck, &
         dx = dx,       dy = dy,            &
         thck_effective = thck_effective,   &
         partial_cf_mask = partial_cf_mask, &
         full_mask = full_mask)

    if (verbose_calving) then
       call point_diag(thck, 'Eigencalving, thck (m)', itest, jtest, rtest, 7, 7, '(f10.3)')
       call point_diag(eps_eigen1*scyr, 'eps1 (1/yr)',    itest, jtest, rtest, 7, 7, '(f10.6)')
       call point_diag(eps_eigen2*scyr, 'eps2 (1/yr)',    itest, jtest, rtest, 7, 7, '(f10.6)')
       call point_diag(tau_eigen1,      'tau1 (Pa)',      itest, jtest, rtest, 7, 7, '(f10.2)')
       call point_diag(tau_eigen2,      'tau2 (Pa)',      itest, jtest, rtest, 7, 7, '(f10.2)')
    endif

    ! Use the eigenvalues to compute a lateral calving rate.

    lateral_rate(:,:) = 0.0d0
    iterate_calving = .true.
    count = 0

    where (calving_front_mask == 1)
       dt_calving = dt
    elsewhere
       dt_calving = 0.0d0
    endwhere

    !TODO - Much of the following code is repeated in damage-based calving.
    !       Put duplicate code in subroutines?

    do while (iterate_calving)

       count = count + 1
       if (verbose_calving .and. this_rank == rtest) then
          write(6,*) 'Iterate calving, count =', count
       endif

       iterate_calving = .false.
       thck_old = thck   ! save current value of thck

       ! For each cell, compute an effective ice thickness for calving.
       ! Away from the calving front, we have thck_effective = thck.
       ! For CF cells, thck_effective is a weighted average of the thickness in a cell and its thickest neighbor.
       ! We require thck_effective >= calving%minthck, as well as thck_effective >= thck.

       if (verbose_calving) then
          call point_diag(thck_effective, 'Compute thck_effective (m)', itest, jtest, rtest, 7, 7)
          call point_diag(thck, 'thck (m)', itest, jtest, rtest, 7, 7)
          call point_diag(partial_cf_mask, 'partial_cf_mask', itest, jtest, rtest, 7, 7)
          call point_diag(full_mask, 'full_mask', itest, jtest, rtest, 7, 7)
          call point_diag(calving_front_mask, 'calving_front_mask', itest, jtest, rtest, 7, 7)
       endif

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (calving_front_mask(i,j) == 1 .and. dt_calving(i,j) > 0.0d0) then

                if (verbose_calving .and. this_rank == rtest .and. count > 1) then
                   write(6,*) 'Next round of calving: count, i, j, dt_calving (yr):', &
                        count, i, j, dt_calving(i,j)/scyr
                endif

                ! Compute a rate factor to account for the length of the calving front.
                ! A typical CF cell will have ice-free ocean on just one edge.
                !  Cells with 2 or 3 ocean neighbors have a longer CF and will calve more quickly.
                ! Also, if a partial_CF cell (H = H1) shares an edge with a thinner partial_CF cell (H = H2 < H1),
                !  we include a term proportional to (1 - H2/H1).
                ! This can be viewed as a subgrid parameterization of calving front length.
                ! It helps to even out the ice thickness along the CF, inhibiting the formation
                !  of notches and protusions.

                ocean_neighbor_count = &
                     ocean_mask(i-1,j) + ocean_mask(i+1,j) + ocean_mask(i,j-1) + ocean_mask(i,j+1)
                lateral_rate_factor = real(ocean_neighbor_count, dp)

                ! loop over edge neighbors, accounting for partial CF cells
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if (ii == i .neqv. jj == j) then   ! edge neighbor (neqv = exclusive or)
                         if (partial_cf_mask(ii,jj) == 1 .and. thck_old(i,j) > thck_old(ii,jj)) then
                            lateral_rate_factor = &
                                 lateral_rate_factor + (1.d0 - thck_old(ii,jj)/thck_old(i,j))
                         endif
                      endif   ! edge neighbor
                   enddo   ! ii
                enddo   ! jj

                !WHL - debug
                if (verbose_calving .and. this_rank == rtest .and. abs(i-itest) <= 1 .and. abs(j-jtest) <= 1) then
                   print*, 'r, i, j, lateral_rate_factor:', rtest, i, j, lateral_rate_factor
                endif

                ! Compute a lateral calving rate (m/s)
                ! Note: tau has the same units as rhoi*g*H, so the eigencalving constants
                !       have units of (1/s).
                !TODO - Add an option to compute the calving rate based on eps_eigen1 and eps_eigen2

                lateral_rate(i,j) = lateral_rate_factor *                                     &
                     (eigenconstant1 * tau_eigen1(i,j) + eigenconstant2 * tau_eigen2(i,j))    &
                     / (rhoi * grav * stress_scale_thck)

                ! Convert the lateral rate to a thinning rate.
                ! For the conversion, assume L * H_eff = dH/dt * l_edge,
                !  where L is the lateral rate and H_eff is the effective thickness of calving ice.
                ! The idea is that the volume removed is equal to the volume that would be removed
                !  by lateral calving at a CF with H = H_eff.

                thinning_rate = lateral_rate(i,j) * thck_effective(i,j) / sqrt(dx*dy)
                dthck = thinning_rate * dt_calving(i,j)  ! m

                ! Compute the new ice thickness

                if (dthck > thck(i,j)) then
                   ! calve the full column and compute the time remaining for more calving
                   dt_calving(i,j) = dt_calving(i,j) * (1.0d0 - thck(i,j)/dthck)
                   iterate_calving = .true.
                   if (verbose_calving .and. this_rank == rtest) then
                      write(6,*) 'Calving time remaining (yr):', i, j, dt_calving(i,j)/scyr
                   endif
                   calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                   thck(i,j) = 0.0d0
                else
                   ! calve part of the column
                   thck(i,j) = thck(i,j) - dthck
                   calving_thck(i,j) = calving_thck(i,j) + dthck
                   dt_calving(i,j) = 0.0d0
                endif

             else   ! no calving

                dt_calving(i,j) = 0.0d0

             endif   ! calving_front cell with damage > threshold

          enddo   ! i
       enddo   ! j

       call parallel_halo(thck, parallel)
       call parallel_halo(calving_thck, parallel)

       if (verbose_calving) then
          call point_diag(lateral_rate*scyr, 'lateral rate (m/yr)', itest, jtest, rtest, 7, 7)
       endif

       if (iterate_calving) then

          ! Recompute the CF mask
          ! This could result in new CF cells that were previously in the interior.
          ! Note: thck_effective keeps the value computed above.

          calving_front_mask_old = calving_front_mask

          call glissade_get_masks(&
               nx,            ny,             &
               parallel,                      &
               thck,          topg,           &
               eus,           thklim,         &
               ice_mask,                      &
               floating_mask = floating_mask, &
               ocean_mask = ocean_mask,       &
               land_mask = land_mask)

          !TODO - Call without thck_effective argument
          call glissade_calving_front_mask(&
               nx,            ny,                 &
               which_ho_calving_front,            &
               parallel,                          &
               thck,          topg,               &
               eus,                               &
               ice_mask,      floating_mask,      &
               ocean_mask,    land_mask,          &
               calving_front_mask,                &
               thck_calving_front)

          ! Save the old dt_calving values and compute new values
          ! Loop through cells; if a new CF cell has ice-free neighbors with dt_calving > 0,
          !  then use this time to calve.  If it has multiple such neighbors, then use the largest value.

          dt_calving_old = dt_calving
          dt_calving = 0.0d0

          iterate_calving = .false.  ! reset to false, in case there are no calving-ready neighbors
          do j = 2, ny-1
             do i = 2, nx-1
                if (calving_front_mask_old(i,j) == 0 .and. calving_front_mask(i,j) == 1) then
                   dt_calving(i,j) = max(dt_calving_old(i-1,j), dt_calving_old(i+1,j), &
                                         dt_calving_old(i,j-1), dt_calving_old(i,j+1))
                   iterate_calving = .true.
                   !WHL - debug
                   if (dt_calving(i,j) > 0.0d0) then
                      write(6,*) 'rank, i, j, new dt_calving (yr):', this_rank, i, j, dt_calving(i,j)/scyr
                   endif
                endif
             enddo
          enddo

       endif   ! iterate_calving

       if (verbose_calving) then
          call point_diag(thck, 'Done iterating, thck', itest, jtest, rtest, 7, 7)
          call point_diag(calving_thck, 'calving_thck', itest, jtest, rtest, 7, 7)
       endif

    enddo   ! do while iterate_calving = T

    if (verbose_calving .and. this_rank == rtest) then
       write(6,*) ' '
       write(6,*) 'Done in eigencalving'
    endif

  end subroutine eigen_calving

!---------------------------------------------------------------------------

  subroutine damage_based_calving(&
       nx,       ny,       nz,                    &
       dx,                 dy,                    &
       sigma,              dt,                    &
       itest,    jtest,    rtest,                 &
       parallel,                                  &
       which_ho_calving_front,                    &
       velnorm_mean,                              &
       thck,               topg,                  &
       eus,                thklim,                &
       tau_eigen1,         tau_eigen2,            &
       eps_eigen1,         eps_eigen2,            &
       damage_constant1,   damage_constant2,      &
       damage_threshold,   calving_minthck,       &
       damage,                                    &
       lateral_rate,       calving_thck)

    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask
    use glissade_grid_operators, only: glissade_unstagger

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny, nz,             & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt                        ! time step (s)

    real(dp), dimension(nz), intent(in) ::&
         sigma                     ! vertical sigma coordinate

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    integer, intent(in) :: &
         which_ho_calving_front    ! = 1 for subgrid calving-front scheme, else = 0

    real(dp), dimension(nx-1,ny-1), intent(in) :: &
         velnorm_mean              ! mean ice speed at vertices (m/s)

    real(dp), dimension(nx,ny), intent(in) :: &
         topg,                   & ! bed topography (m)
         tau_eigen1, tau_eigen2, & ! eigenvalues of the horizontal stress tensor (Pa)
         eps_eigen1, eps_eigen2    ! eigenvalues of the horizontal strain rate tensor (1/s)

    real(dp), intent(in) ::  &
         eus,                    & ! eustatic sea level (m)
         thklim                    ! minimum thickness for dynamically active grounded ice (m)


    real(dp), intent(in) :: &
         damage_constant1,       & ! rate of change of damage (1/s) proportional to tau1
         damage_constant2,       & ! rate of change of damage (1/s) proportional to tau2
         damage_threshold,       & ! calving begins for damage > damage_threshold
         calving_minthck           ! min effective thickness (m) for CF cells

    real(dp), dimension(nz-1,nx,ny), intent(inout) :: &
         damage                    ! 3D damage tracer, 0 > damage < 1

    real(dp), dimension(nx,ny), intent(out) :: &
         lateral_rate              ! lateral rate of calving (m/s)

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m)
         calving_thck              ! thickness lost due to calving (m)

    ! local variables

    integer :: i, j, k, ii, jj

    integer, dimension(nx,ny)  ::  &
         ice_mask,               & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,          & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,             & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,              & ! = 1 where topg is at or above sea level, else = 0
         calving_front_mask,     & ! = 1 where ice is floating and borders at least one ocean cell, else = 0
         calving_front_mask_old

    integer, dimension(nx,ny) :: &
         partial_cf_mask,        & ! = 1 for partially filled CF cells (thck < thck_effective), else = 0
         full_mask                 ! = 1 for ice-filled cells that are not partial_cf cells, else = 0

    integer, dimension(nx-1,ny-1) :: &
         vmask                     ! = 1 for vertices with ice speed (i.e., vertices of active cells)

    real(dp), dimension(nx,ny) :: &
         speed,                  & ! 2D ice speed averaged to cell centers (m/s)
         damage_column,          & ! 2D vertically integrated scalar damage parameter
         thck_calving_front,     & ! effective ice thickness at the calving front (m)
         thck_effective,         & ! effective thickness for calving (m)
         thck_old,               & ! old value of thck (m)
         dt_calving,             & ! time remaining for calving (s)
         dt_calving_old            ! old value of dt_calving

    real(dp) :: &
         stress_scale,           & ! stress scale in prognostic damage equation (Pa)
         d_damage_dt,            & ! rate of change of damage scalar (1/s)
         damage_frac,            & ! fraction of maximum thinning rate, based on damage
         ocean_neighbor_count,   & ! number of edge neighbors that are ice-free ocean
         lateral_rate_factor,    & ! multiplying factor to account for a CF longer than one edge
         thinning_rate,          & ! vertical thinning rate (m/s)
         dthck                     ! thickness change (m)

    logical, parameter :: &
         damage_scale_stress_by_thickness = .true.   ! if true, stress_scale is proportional to H
                                                     ! if false, stress_scale is independent of local H
    real(dp), parameter :: &
         stress_scale_thck = 100.d0     ! ice thickness (m) in stress_scale = rhoi*g*H

    logical, parameter :: &
         stress_based_damage = .true.   ! true if dD/dt is based on tau1 and tau2, false if based on eps1 and eps2

    logical :: iterate_calving          ! Iterate the calving until dt_calving = 0 for all cells
    integer :: count

    ! Not sure these parallel updates are needed here
    call parallel_halo(tau_eigen1, parallel)
    call parallel_halo(tau_eigen2, parallel)
    call parallel_halo(eps_eigen1, parallel)
    call parallel_halo(eps_eigen2, parallel)

    ! Compute masks needed for calving
    ! Need a calving_front_mask; calving/thinning is applied only to cells at the calving front.
    ! For each CF cell, also compute an effective ice thickness for calving.

    call glissade_get_masks(&
         nx,            ny,             &
         parallel,                      &
         thck,          topg,           &
         eus,           thklim,         &
         ice_mask,                      &
         floating_mask = floating_mask, &
         ocean_mask = ocean_mask,       &
         land_mask = land_mask)

    call glissade_calving_front_mask(&
         nx,            ny,                 &
         which_ho_calving_front,            &
         parallel,                          &
         thck,          topg,               &
         eus,                               &
         ice_mask,      floating_mask,      &
         ocean_mask,    land_mask,          &
         calving_front_mask,                &
         thck_calving_front,                &
         calving_minthck = calving_minthck, &
         dx = dx,       dy = dy,            &
         thck_effective = thck_effective,   &
         partial_cf_mask = partial_cf_mask, &
         full_mask = full_mask)

    ! Compute the ice speed at cell centers, averaged from neighboring vertices.
    ! Include in the average only vertices with nonzero speeds (i.e., ice present)
    where (velnorm_mean > 0.0d0)
       vmask = 1
    elsewhere
       vmask = 0
    endwhere

    call glissade_unstagger(&
         nx,            ny,    &
         velnorm_mean,  speed, &
         vmask,         stagger_margin_in = 1)

    ! Compute the vertically integrated damage in each column (diagnostic only)
    damage_column(:,:) = 0.0d0
    do j = 1, ny
       do i = 1, nx
          do k = 1, nz-1
             damage_column(i,j) = damage_column(i,j) + damage(k,i,j) * (sigma(k+1) - sigma(k))
          enddo
       enddo
    enddo

    if (verbose_calving) then
       call point_diag(speed*scyr, 'Damage-based calving, speed (m/yr)', &
            itest, jtest, rtest, 7, 7, '(f10.2)')
       call point_diag(eps_eigen1*scyr, 'eps1 (1/yr)',    itest, jtest, rtest, 7, 7, '(f10.6)')
       call point_diag(eps_eigen2*scyr, 'eps2 (1/yr)',    itest, jtest, rtest, 7, 7, '(f10.6)')
       call point_diag(tau_eigen1,      'tau1 (Pa)',      itest, jtest, rtest, 7, 7, '(f10.2)')
       call point_diag(tau_eigen2,      'tau2 (Pa)',      itest, jtest, rtest, 7, 7, '(f10.2)')
       call point_diag(damage_column,   'initial damage', itest, jtest, rtest, 7, 7, '(f10.6)')
    endif

    ! Prognose changes in damage.
    ! For now, this is done using a simple scheme based on the horizontal stress tensor.
    ! At some point, we may want to prognose damage in a way that depends on other factors such as mass balance.
    ! Note: Damage is formally a 3D field, which makes it easier to advect, even though
    !       (in the current scheme) the damage source term is uniform in each column.

    ! Compute the change in damage in each cell
    do j = 2, ny-1
       do i = 1, nx-1
          if (floating_mask(i,j) == 1) then
             if (stress_based_damage) then

                ! Two cases:
                ! (1) d_damage_dt is proportional to tau_eff/H
                ! (2) d_damage_dt is proportional to tau_eff
                ! This damage_constant has units of s^{-1}

                if (damage_scale_stress_by_thickness) then
                   stress_scale = rhoi*grav * max(thck(i,j), stress_scale_thck)
                else
                   stress_scale = rhoi*grav * stress_scale_thck
                endif

                d_damage_dt = (damage_constant1 * max(tau_eigen1(i,j), 0.0d0) +  &
                               damage_constant2 * max(tau_eigen2(i,j), 0.0d0))   &
                               / stress_scale

             else   ! strain-based damage

                ! Note: The damage constant is a scalar; eps and d_damage_dt have units of 1/s
                d_damage_dt = (damage_constant1*scyr * max(eps_eigen1(i,j), 0.0d0) +  &
                               damage_constant2*scyr * max(eps_eigen2(i,j), 0.0d0))

             endif  ! stress or strain-based

             damage(:,i,j) = damage(:,i,j) + d_damage_dt * dt
             damage(:,i,j) = min(damage(:,i,j), 1.0d0)
             damage(:,i,j) = max(damage(:,i,j), 0.0d0)

          else    ! set damage to zero elsewhere (grounded or ice-free)
             damage(:,i,j) = 0.0d0
          endif   ! floating_mask

       enddo   ! i
    enddo   ! j

    ! Compute the vertically integrated damage in each column.
    damage_column(:,:) = 0.0d0
    do j = 1, ny
       do i = 1, nx
          do k = 1, nz-1
             damage_column(i,j) = damage_column(i,j) + damage(k,i,j) * (sigma(k+1) - sigma(k))
          enddo
       enddo
    enddo

    if (verbose_calving) then
       call point_diag(damage_column, 'new damage', itest, jtest, rtest, 7, 7, '(f10.6)')
    endif

    ! Main damage-calving loops
    ! Loop over locally owned cells
    ! Calving occurs only in CF cells: cells with one or more ocean neighbors.
    ! If time remains after one or more columns have calved completely, iterate until dt_calving = 0 everywhere.

    lateral_rate(:,:) = 0.0d0
    iterate_calving = .true.
    count = 0

    where (calving_front_mask == 1)
       dt_calving = dt
    elsewhere
       dt_calving = 0.0d0
    endwhere

    do while (iterate_calving)

       count = count + 1
       if (verbose_calving .and. this_rank == rtest) then
          write(6,*) 'Iterate calving, count =', count
       endif

       iterate_calving = .false.
       thck_old = thck   ! save current value of thck

       if (verbose_calving) then
          call point_diag(thck_effective, 'Compute thck_effective (m)', itest, jtest, rtest, 7, 7)
          call point_diag(thck, 'thck (m)', itest, jtest, rtest, 7, 7)
          call point_diag(partial_cf_mask, 'partial_cf_mask', itest, jtest, rtest, 7, 7)
          call point_diag(full_mask, 'full_mask', itest, jtest, rtest, 7, 7)
          call point_diag(calving_front_mask, 'calving_front_mask', itest, jtest, rtest, 7, 7)
       endif

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (calving_front_mask(i,j) == 1 .and. damage_column(i,j) > damage_threshold  &
                  .and. dt_calving(i,j) > 0.0d0) then

                if (verbose_calving .and. this_rank == rtest .and. count > 1) then
                   write(6,*) 'Next round of calving: count, i, j, dt_calving (yr):', &
                        count, i, j, dt_calving(i,j)/scyr
                endif

                ! Compute a damage fraction, increasing linearly from 0 (at a damage threshold) to 1.
                damage_frac = (damage_column(i,j) - damage_threshold) / (1.0d0 - damage_threshold)

                ! Compute a rate factor to account for the length of the calving front.
                ! A typical CF cell will have ice-free ocean on just one edge.
                !  Cells with 2 or 3 ocean neighbors have a longer CF and will calve more quickly.
                ! Also, if a partial_CF cell (H = H1) shares an edge with a thinner partial_CF cell (H = H2 < H1),
                !  we include a term proportional to (1 - H2/H1).
                ! This can be viewed as a subgrid parameterization of calving front length.
                ! It helps to even out the ice thickness along the CF, inhibiting the formation
                !  of notches and protusions.

                ocean_neighbor_count = &
                     ocean_mask(i-1,j) + ocean_mask(i+1,j) + ocean_mask(i,j-1) + ocean_mask(i,j+1)
                lateral_rate_factor = real(ocean_neighbor_count, dp)

                ! loop over edge neighbors, accounting for partial CF cells
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if (ii == i .neqv. jj == j) then   ! edge neighbor (neqv = exclusive or)
                         if (partial_cf_mask(ii,jj) == 1 .and. thck_old(i,j) > thck_old(ii,jj)) then
                            lateral_rate_factor = &
                                 lateral_rate_factor + (1.d0 - thck_old(ii,jj)/thck_old(i,j))
                         endif
                      endif   ! edge neighbor
                   enddo   ! ii
                enddo   ! jj

                if (verbose_calving .and. this_rank == rtest .and. abs(i-itest) <= 3 .and. abs(j-jtest) <= 3) then
                   print*, 'r, i, j, lateral_rate_factor:', rtest, i, j, lateral_rate_factor
                endif

                ! Compute a lateral calving rate (m/s)
                ! The line commented out assumes a prescribed maximum rate instead of the local speed.
                ! Using the local speed will generally halt CF advance when damage_frac ~ 1.
                lateral_rate(i,j) = lateral_rate_factor * damage_frac * speed(i,j)

                ! Convert the lateral rate to a thinning rate.
                ! For the conversion, assume L * H_eff = dH/dt * l_edge,
                !  where L is the lateral rate and H_eff is the effective thickness of calving ice.
                ! The idea is that the volume removed is equal to the volume that would be removed
                !  by lateral calving at a CF with H = H_eff.

                thinning_rate = lateral_rate(i,j) * thck_effective(i,j) / sqrt(dx*dy)
                dthck = thinning_rate * dt_calving(i,j)  ! m

                ! Compute the new ice thickness

                if (dthck > thck(i,j)) then
                   ! calve the full column and compute the time remaining for more calving
                   dt_calving(i,j) = dt_calving(i,j) * (1.0d0 - thck(i,j)/dthck)
                   iterate_calving = .true.
                   if (verbose_calving .and. this_rank == rtest) then
                      write(6,*) 'Calving time remaining (yr):', i, j, dt_calving(i,j)/scyr
                   endif
                   calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                   thck(i,j) = 0.0d0
                else
                   ! calve part of the column
                   thck(i,j) = thck(i,j) - dthck
                   calving_thck(i,j) = calving_thck(i,j) + dthck
                   dt_calving(i,j) = 0.0d0
                endif

             else   ! no calving

                dt_calving(i,j) = 0.0d0

             endif   ! calving_front cell with damage > threshold

          enddo   ! i
       enddo   ! j

       call parallel_halo(thck, parallel)
       call parallel_halo(calving_thck, parallel)

       if (verbose_calving) then
          call point_diag(lateral_rate*scyr, 'lateral rate (m/yr)', itest, jtest, rtest, 7, 7)
       endif

       if (iterate_calving) then

          ! Recompute the CF mask
          ! This could result in new CF cells that were previously in the interior.
          ! Note: thck_effective keeps the value computed above.

          calving_front_mask_old = calving_front_mask

          call glissade_get_masks(&
               nx,            ny,             &
               parallel,                      &
               thck,          topg,           &
               eus,           thklim,         &
               ice_mask,                      &
               floating_mask = floating_mask, &
               ocean_mask = ocean_mask,       &
               land_mask = land_mask)

          !TODO - Call without thck_calving_front argument
          call glissade_calving_front_mask(&
               nx,            ny,                 &
               which_ho_calving_front,            &
               parallel,                          &
               thck,          topg,               &
               eus,                               &
               ice_mask,      floating_mask,      &
               ocean_mask,    land_mask,          &
               calving_front_mask,                &
               thck_calving_front)

          ! Save the old dt_calving values and compute new values
          ! Loop through cells; if a new CF cell has ice-free neighbors with dt_calving > 0,
          !  then use this time to calve.  If it has multiple such neighbors, then use the largest value.

          dt_calving_old = dt_calving
          dt_calving = 0.0d0

          iterate_calving = .false.  ! reset to false, in case there are no calving-ready neighbors
          do j = 2, ny-1
             do i = 2, nx-1
                if (calving_front_mask_old(i,j) == 0 .and. calving_front_mask(i,j) == 1  &
                     .and. damage_column(i,j) > damage_threshold) then   ! new CF cell that can calve
                   dt_calving(i,j) = max(dt_calving_old(i-1,j), dt_calving_old(i+1,j), &
                                         dt_calving_old(i,j-1), dt_calving_old(i,j+1))
                   iterate_calving = .true.
                   !WHL - debug
                   if (dt_calving(i,j) > 0.0d0) then
                      write(6,*) 'rank, i, j, new dt_calving (yr):', this_rank, i, j, dt_calving(i,j)/scyr
                   endif
                endif
             enddo
          enddo

       endif   ! iterate_calving

       if (verbose_calving) then
          call point_diag(thck, 'Done iterating, thck', itest, jtest, rtest, 7, 7)
          call point_diag(calving_thck, 'calving_thck', itest, jtest, rtest, 7, 7)
       endif

    enddo   ! do while iterate_calving = T

    if (verbose_calving .and. this_rank == rtest) then
       write(6,*) ' '
       write(6,*) 'Done in damage-based calving'
    endif

  end subroutine damage_based_calving

!---------------------------------------------------------------------------

  subroutine glissade_cull_calving_front(&
       nx,           ny,          &
       parallel,                  &
       itest, jtest, rtest,       &
       thck,         topg,        &
       eus,          thklim,      &
       which_ho_calving_front,    &
       ncull_calving_front,       &
       calving_thck)

    ! Optionally, remove all cells currently at the calving front are removed.
    ! This subroutine would typically be called at initialization, if at all.
    ! Culling can removed long, skinny floating peninsulas that can be dynamically unstable.
    !  Without this step, peninsulas up to two cells thick (with calving-front cells on each side)
    !  will typically be removed as icebergs (because there is no path back to grounded ice through active cells).
    ! With one round of culling, peninsulas up to four cells thick will be removed (two outer layers
    !  during the preliminary step, followed by two inner layers on the remove_iceberg step).
    ! If necessary, culling can be repeated to remove peninsulas with a thickness of 6 layers, 8 layers, etc.

    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask

    integer, intent(in) :: nx, ny                       !> horizontal grid dimensions
    type(parallel_type), intent(in) :: parallel         !> info for parallel communication
    integer, intent(in) :: itest, jtest, rtest          !> coordinates of diagnostic point

    real(dp), dimension(nx,ny), intent(inout) :: thck   !> ice thickness
    real(dp), dimension(nx,ny), intent(in)    :: topg   !> present bedrock topography
    real(dp), intent(in)    :: eus                      !> eustatic sea level
    real(dp), intent(in)    :: thklim                   !> minimum thickness for dynamically active grounded ice
    integer, intent(in)     :: which_ho_calving_front   !> = 1 for subgrid calving-front scheme, else = 0
    integer, intent(in) :: &
         ncull_calving_front           !> number of times to cull calving_front cells at initialization

    real(dp), dimension(nx,ny), intent(inout) :: calving_thck   !> thickness lost due to calving in each grid cell;
                                                              !> on output, includes ice that is culled here
    ! local variables

    integer :: i, j, n

    integer,  dimension(nx,ny) ::  &
         ice_mask,           & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,      & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,         & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,          & ! = 1 where topg is at or above sea level, else = 0
         calving_front_mask    ! = 1 where ice is floating and borders the ocean, else = 0

    real(dp),  dimension(nx,ny) ::  &
         thck_calving_front    ! effective ice thickness at the calving front

    do n = 1, ncull_calving_front

       ! calculate masks
       ! Note: Passing in thklim = 0.0 does not work because it erroneously counts thin floating cells as active.
       !       Then the algorithm can fail to identify floating regions that should be removed
       !       (since they are separated from any active cells).

       call glissade_get_masks(nx,            ny,             &
                               parallel,                      &
                               thck,          topg,           &
                               eus,           thklim,         &
                               ice_mask,                      &
                               floating_mask = floating_mask, &
                               ocean_mask = ocean_mask,       &
                               land_mask = land_mask)

       call glissade_calving_front_mask(nx,            ny,                 &
                                        which_ho_calving_front,            &
                                        parallel,                          &
                                        thck,          topg,               &
                                        eus,                               &
                                        ice_mask,      floating_mask,      &
                                        ocean_mask,    land_mask,          &
                                        calving_front_mask,                &
                                        thck_calving_front)

       if (main_task) then
          call write_log ('cull_calving_front: Removing ice from calving_front cells')
          write(6,*) 'cull_calving_front: Removing ice from calving_front cells'
       endif

       if (verbose_calving) then
          call point_diag(calving_front_mask, 'calving_front_mask for culling', &
               itest, jtest, rtest, 7, 7)
          call point_diag(thck, 'thck (before CF culling)', itest, jtest, rtest, 7, 7)
       endif

       do j = 1, ny
          do i = 1, nx
             if (calving_front_mask(i,j) == 1) then
                calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                thck(i,j) = 0.0d0
             endif
          enddo
       enddo

       if (verbose_calving) then
          call point_diag(thck, 'thck (after CF culling)', itest, jtest, rtest, 7, 7)
       endif

    enddo  ! ncull_calving_front

  end subroutine glissade_cull_calving_front

!---------------------------------------------------------------------------

  subroutine glissade_remove_icebergs(&
       nx,           ny,            &
       parallel,                    &
       itest, jtest, rtest,         &
       thck,                        &
       f_ground_cell,               &
       ice_mask,                    &
       land_mask,                   &
       active_ice_mask,             &
       calving_thck)

    ! Remove any icebergs. 
        
    ! The algorithm is as follows:
    ! (1) Mark all cells with ice (either active or inactive) with the initial color.
    !     Mark other cells with the boundary color.
    ! (2) Seed the fill by giving active grounded cells the fill color.
    ! (3) Recursively fill all cells that are connected to filled cells by a path
    !     that passes through active cells only.
    ! (4) Repeat the recursion as necessary to spread the fill to adjacent processors.
    ! (5) Once the fill is done, any floating cells that still have the initial color
    !     are considered to be icebergs and are removed.
    !
    ! Notes:
    ! (1) Grounded cells must have f_ground_cell > f_ground_threshold to seed the fill.
    ! (2) The recursive fill applies to edge neighbors, not corner neighbors.
    !     The path back to grounded ice must go through edges, not corners.
    ! (3) Inactive cells can be filled (if adjacent to active cells), but
    !     do not further spread the fill.
    ! (4) Should have thklim > 0.  With a limit of 0.0, very thin floating cells
    !     can be wrongly counted as active, and icebergs can be missed.
    ! (5) Land-based cells that still have the initial color are not marked as icebergs.

    use glissade_masks, only: glissade_fill_with_buffer, initial_color, fill_color, boundary_color

    integer, intent(in) :: nx, ny                                !> horizontal grid dimensions
    type(parallel_type), intent(in) :: parallel                  !> info for parallel communication
    integer, intent(in) :: itest, jtest, rtest                   !> coordinates of diagnostic point

    real(dp), dimension(nx,ny), intent(inout) :: thck            !> ice thickness
    real(dp), dimension(nx,ny), intent(in)    :: f_ground_cell   !> grounded fraction in each grid cell
    integer,  dimension(nx,ny), intent(in)    :: ice_mask        !> = 1 where ice is present (thck > thklim), else = 0
    integer,  dimension(nx,ny), intent(in)    :: land_mask       !> = 1 where topg - eus >= 0, else = 0
    integer,  dimension(nx,ny), intent(in)    :: active_ice_mask !> = 1 for dynamically active cells
    real(dp), dimension(nx,ny), intent(inout) :: calving_thck    !> thickness lost due to calving in each grid cell;
                                                                 !> on output, includes ice in icebergs
    ! local variables

    integer :: i, j, iter

    integer :: &
         max_iter,             & ! max(ewtasks, nstasks)
         local_count,          & ! local counter for filled values
         global_count,         & ! global counter for filled values
         global_count_save       ! globalcounter for filled values from previous iteration

    integer,  dimension(nx,ny) ::  &
         color                 ! integer 'color' for identifying icebergs

    real(dp),  dimension(nx,ny) ::  &
         thck_calving_front    ! effective ice thickness at the calving front

    !TODO - Make this a config parameter?
    real(dp), parameter :: &   ! threshold for counting cells as grounded
         f_ground_threshold = 0.10d0

    if (verbose_calving) then
       call point_diag(thck, 'Remove icebergs, thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(active_ice_mask, 'active_ice_mask', itest, jtest, rtest, 7, 7)
       call point_diag(f_ground_cell, 'f_ground_cell', itest, jtest, rtest, 7, 7)
    endif
    
    ! Initialize iceberg removal
    ! Note: Any cell with ice, active or inactive, receives the initial color.
    !       Inactive cells can later receive the fill color (if adjacent to active cells)
    !        but cannot further spread the fill color.
    !       This protects inactive calving-front cells from removal, as desired.

    do j = 1, ny
       do i = 1, nx
          if (thck(i,j) > 0.0d0) then
             color(i,j) = initial_color
          else
             color(i,j) = boundary_color
          endif
       enddo
    enddo

    ! Loop through cells, identifying active cells with grounded ice.
    ! Fill each grounded cell and then recursively fill active neighbor cells, whether grounded or not.
    ! We may have to do this several times to incorporate connections between neighboring processors.

    max_iter = max(parallel%ewtasks, parallel%nstasks)
    global_count_save = 0

    do iter = 1, max_iter

       if (iter == 1) then   ! identify grounded cells that can seed the fill

          do j = 1, ny
             do i = 1, nx

                ! Active cells that are firmly grounded can seed the fill.
                ! It is somewhat arbitrary how to define "firmly grounded", but here we require
                !  that f_ground_cell exceeds a threshold value defined above.
                ! Note: If running without a GLP, then f_ground_cell is binary, either 0 or 1.

                if (active_ice_mask(i,j) == 1 .and. f_ground_cell(i,j) > f_ground_threshold) then  ! grounded ice

                   if (color(i,j) /= boundary_color .and. color(i,j) /= fill_color) then

                      ! assign the fill color to this cell, and recursively fill neighbor cells
                      call glissade_fill_with_buffer(nx,    ny,    &
                                                     i,     j,     &
                                                     color, active_ice_mask)

                   endif

                endif
             enddo
          enddo

       else  ! count > 1

          ! Check for halo cells that were just filled on neighbor processors
          ! Note: In order for a halo cell to seed the fill on this processor, it must not only have the fill color,
          !       but also must be an active cell.

          call parallel_halo(color, parallel)

          ! west halo layer
          i = nhalo
          do j = 1, ny
             if (color(i,j) == fill_color .and. active_ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(nx,    ny,    &
                                               i+1,   j,     &
                                               color, active_ice_mask)
             endif
          enddo

          ! east halo layers
          i = nx - nhalo + 1
          do j = 1, ny
             if (color(i,j) == fill_color .and. active_ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(nx,    ny,    &
                                               i-1,   j,     &
                                               color, active_ice_mask)
             endif
          enddo

          ! south halo layer
          j = nhalo
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. active_ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(nx,    ny,    &
                                               i,     j+1,   &
                                               color, active_ice_mask)
             endif
          enddo

          ! north halo layer
          j = ny-nhalo+1
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. active_ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(nx,    ny,    &
                                               i,     j-1,   &
                                               color, active_ice_mask)
             endif
          enddo

       endif  ! count = 1

       local_count = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (color(i,j) == fill_color) local_count = local_count + 1
          enddo
       enddo

       !WHL - If running a large problem, may want to reduce the frequency of this global sum
       global_count = parallel_reduce_sum(local_count)

       if (global_count == global_count_save) then
          if (verbose_calving .and. main_task) &
               write(6,*) 'Fill converged: iter, global_count =', iter, global_count
          exit
       else
          if (verbose_calving .and. main_task) &
               write(6,*) 'Convergence check: iter, global_count =', iter, global_count
          global_count_save = global_count
       endif

    enddo  ! count

    ! Icebergs are cells that still have the initial color and are not on land.
    ! Remove ice in these cells, adding it to the calving field.
    ! Note: There is an exception for cells that are
    !       (1) adjacent to at least one ice-covered cell (sharing an edge), and
    !       (2) connected diagonally to an active cell with the fill color.
    !       Such cells are considered part of the inactive calving front and are
    !        allowed to continue filling instead of calving.

    do j = 2, ny-1
       do i = 2, nx-1
          if (color(i,j) == initial_color .and. land_mask(i,j) == 0) then
             if (  ( color(i-1,j+1)==fill_color .and. active_ice_mask(i-1,j+1)==1 .and. &
                       (ice_mask(i-1,j)==1 .or. ice_mask(i,j+1)==1) ) &
              .or. ( color(i+1,j+1)==fill_color .and. active_ice_mask(i+1,j+1)==1 .and. &
                       (ice_mask(i+1,j)==1 .or. ice_mask(i,j+1)==1) ) &
              .or. ( color(i-1,j-1)==fill_color .and. active_ice_mask(i-1,j-1)==1 .and. &
                       (ice_mask(i-1,j)==1 .or. ice_mask(i,j-1)==1) ) &
              .or. ( color(i+1,j-1)==fill_color .and. active_ice_mask(i+1,j-1)==1 .and. &
                       (ice_mask(i+1,j)==1 .or. ice_mask(i,j-1)==1) ) ) then
                ! do nothing; this cell is part of the inactive calving front
             else  ! not part of the inactive calving front; calve as an iceberg
                calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                thck(i,j) = 0.0d0
                !TODO - Also handle tracers?  E.g., set damage(:,i,j) = 0.d0?
             endif  ! diagonally connected or not
          endif
       enddo
    enddo

    if (verbose_calving) then
       call point_diag(thck, 'After iceberg removal, thck', itest, jtest, rtest, 7, 7)
    endif

  end subroutine glissade_remove_icebergs

!---------------------------------------------------------------------------

  subroutine glissade_remove_isthmuses(&
       nx,           ny,            &
       itest, jtest, rtest,         &
       thck,                        &
       f_ground_cell,               &
       floating_mask,               &
       ocean_mask,                  &
       calving_thck)

    ! Remove any ice isthmuses.
    ! An isthmus is defined as a floating or weakly grounded grid cell with ice-free ocean
    !  or thin floating ice (H < 10 m) on both sides: i.e., in cells (i-1,j) and (i+1,j),
    !  or (i,j-1) and (i,j+1).
    ! When using an ice_fraction_retreat_mask derived from a model (e.g., CESM),
    !  it may be necessary to remove isthmuses to prevent unstable ice configurations,
    !  e.g. a shelf split into two parts connected by a bridge one cell wide.
    ! Isthmus removal should always be followed by iceberg removal.

    integer :: nx, ny                                   !> horizontal grid dimensions
    integer, intent(in) :: itest, jtest, rtest          !> coordinates of diagnostic point

    real(dp), dimension(nx,ny), intent(inout) :: thck            !> ice thickness (m)
    real(dp), dimension(nx,ny), intent(in)    :: f_ground_cell   !> grounded fraction in each grid cell
    integer,  dimension(nx,ny), intent(in)    :: floating_mask   !> = 1 ice is present and floating, else = 0
    integer,  dimension(nx,ny), intent(in)    :: ocean_mask      !> = 1 where topg is below sea level and ice is absent, else = 0
    real(dp), dimension(nx,ny), intent(inout) :: calving_thck    !> thickness (m) lost due to calving in each grid cell;
                                                                 !> on output, includes ice removed from isthmuses

    ! local variables

    integer :: i, j

    integer, dimension(nx,ny) :: &
         ocean_plus_thin_ice_mask         ! = 1 for ocean cells and cells with thin floating ice

    ! Both floating and weakly grounded cells can be identified as isthmuses and removed;
    !  isthmus_f_ground_threshold is used to identify weakly grounded cells.
    real(dp), parameter :: &   ! threshold for counting cells as grounded
         isthmus_f_ground_threshold = 0.50d0

    ! An isthmus cell has ice-free ocean or thin floating ice on each side:
    !  isthmus_f_ground_threshold is used to identify thin floating ice.
    real(dp), parameter :: &   ! threshold (m) for counting floating ice as thin
         isthmus_thck_threshold = 10.0d0

    if (verbose_calving) then
       call point_diag(thck, 'Remove isthmuses, thck', itest, jtest, rtest, 7, 7)
    endif

    ocean_plus_thin_ice_mask = ocean_mask
    where (floating_mask == 1 .and. thck < isthmus_thck_threshold)
       ocean_plus_thin_ice_mask = 1
    endwhere

    do j = 2, ny-1
       do i = 2, nx-1
          if (floating_mask(i,j) == 1 .or. f_ground_cell(i,j) < isthmus_f_ground_threshold) then
             if ( (ocean_plus_thin_ice_mask(i-1,j) == 1 .and. ocean_plus_thin_ice_mask(i+1,j) == 1) .or. &
                  (ocean_plus_thin_ice_mask(i,j-1) == 1 .and. ocean_plus_thin_ice_mask(i,j+1) == 1) ) then
                calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                thck(i,j) = 0.0d0
             endif
          endif
       enddo
    enddo

    if (verbose_calving) then
       call point_diag(thck, 'After isthmus removal, thck', itest, jtest, rtest, 7, 7)
    endif

  end subroutine glissade_remove_isthmuses

!---------------------------------------------------------------------------

  subroutine glissade_limit_cliffs(&
       nx,             ny,              &
       parallel,                        &
       itest,  jtest,  rtest,           &
       dt,                              &
       which_ho_calving_front,          &
       taumax_cliff,   cliff_timescale, &
       thck,           topg,            &
       eus,            thklim,          &
       calving_thck)

    ! Impose a thickness limit on marine ice cliffs.
    ! These are defined as grounded marine-based cells adjacent to inactive calving_front cells or ice-free ocean.
    ! Ice removed from cliffs is added to the calving flux.

    use glissade_masks

    integer, intent(in)  :: nx, ny                      !> horizontal grid dimensions
    type(parallel_type), intent(in) :: parallel         !> info for parallel communication
    integer, intent(in)  :: itest, jtest, rtest         !> coordinates of diagnostic point
    real(dp), intent(in) :: dt                          !> model timestep (s)

    integer, intent(in)     :: which_ho_calving_front   !> = 1 for subgrid calving-front scheme, else = 0
    real(dp), intent(in)    :: taumax_cliff             !> yield stress (Pa) for marine-based ice cliffs
    real(dp), intent(in)    :: cliff_timescale          !> timescale (s) for limiting marine cliff thickness
    real(dp), dimension(nx,ny), intent(inout) :: thck   !> ice thickness (m)
    real(dp), dimension(nx,ny), intent(in)    :: topg   !> present bedrock topography (m)
    real(dp), intent(in)    :: eus                      !> eustatic sea level (m)
    real(dp), intent(in)    :: thklim                   !> minimum thickness for dynamically active grounded ice (m)

    real(dp), dimension(nx,ny), intent(inout) :: &
         calving_thck             !> thickness (m) lost due to calving; on output, includes ice calved at marine cliffs

    ! local variables

    integer :: i, j

    integer, dimension(nx,ny) ::  &
         ice_mask,           & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,      & ! = 1 where ice is present (thck > thklim) and floating, else = 0
         ocean_mask,         & ! = 1 where topg is below sea level and ice is absent, else = 0
         land_mask,          & ! = 1 where topg is at or above sea level, else = 0
         calving_front_mask, & ! = 1 where ice is floating and borders the ocean, else = 0
         active_ice_mask,    & ! = 1 for dynamically active cells
         marine_cliff_mask     ! = 1 where ice is grounded and marine-based and borders at least one
                               !     ocean or inactive calving_front cell, else = 0

    real(dp), dimension(nx,ny) ::  &
         thckmax_cliff,      & ! max stable ice thickness in marine_cliff cells
         thck_calving_front    ! effective ice thickness at the calving front

    real(dp) :: &
         thinning_rate,        & ! vertical thinning rate (m/s)
         dthck,                & ! thickness change (m)
         factor

    ! Update masks, including the marine_cliff mask.
    ! Note: We do not use calving_front_mask or thck_calving_front directly.
    !       But to identify marine cliffs, we use active_ice_mask, which depends on whether there is a subgrid calving front.

    call glissade_get_masks(nx,            ny,             &
                            parallel,                      &
                            thck,          topg,           &
                            eus,           thklim,         &
                            ice_mask,                      &
                            floating_mask = floating_mask, &
                            ocean_mask = ocean_mask,       &
                            land_mask = land_mask)

    call glissade_calving_front_mask(nx,            ny,                 &
                                     which_ho_calving_front,            &
                                     parallel,                          &
                                     thck,          topg,               &
                                     eus,                               &
                                     ice_mask,      floating_mask,      &
                                     ocean_mask,    land_mask,          &
                                     calving_front_mask,                &
                                     thck_calving_front,                &
                                     active_ice_mask = active_ice_mask)

    call glissade_marine_cliff_mask(nx,            ny,                &
                                    ice_mask,      floating_mask,     &
                                    land_mask,     active_ice_mask,   &
                                    marine_cliff_mask)

    call parallel_halo(marine_cliff_mask, parallel)

    if (verbose_calving) then
       call point_diag(marine_cliff_mask, 'Cliff limiting, marine_cliff_mask', itest, jtest, rtest, 7, 7)
       call point_diag(thck, 'thck (m) before limiting', itest, jtest, rtest, 7, 7)
    endif

    thckmax_cliff(:,:) = 0.0d0
    do j = 2, ny-1
       do i = 1, nx-1
          if (marine_cliff_mask(i,j) == 1) then

             ! Compute the max stable ice thickness in the cliff cell.
             ! This is eq. 2.10 in Bassis & Walker (2012)
             factor = taumax_cliff / (rhoi*grav)   ! units are Pa for taumax, m for factor
             thckmax_cliff(i,j) = factor + sqrt(factor**2 + (rhoo/rhoi)*(topg(i,j))**2)  ! m

             ! If thicker than the max stable thickness, then remove some ice and add it to the calving field
             ! Note: By default, cliff_timescale = 0, which means thck is reset to thckmax_cliff each timestep.
             !       Might want to try other values when looking at marine ice cliff instability.
             if (thck(i,j) > thckmax_cliff(i,j)) then

                if (cliff_timescale > 0.0d0) then
                   thinning_rate = (thck(i,j) - thckmax_cliff(i,j)) / cliff_timescale
                   dthck = min(thck(i,j) - thckmax_cliff(i,j), thinning_rate*dt)
                else
                   dthck = thck(i,j) - thckmax_cliff(i,j)
                endif

                thck(i,j) = thck(i,j) - dthck
                calving_thck(i,j) = calving_thck(i,j) + dthck

             endif  ! thck > thckmax_cliff

          endif  ! marine_cliff cell
       enddo   ! i
    enddo   ! j

    if (verbose_calving) then
       call point_diag(thck, 'thck (m) after limiting', itest, jtest, rtest, 7, 7)
       call point_diag(calving_thck, 'calving_thck (m)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine glissade_limit_cliffs

!---------------------------------------------------------------------------

end module glissade_calving

!---------------------------------------------------------------------------
