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
       parallel_halo, parallel_globalindex, &
       parallel_reduce_sum, parallel_reduce_max, parallel_reduce_log_or

  use glimmer_paramets, only: eps08, thk0
  use glimmer_physcon, only: rhoi, rhoo, grav, scyr
  use glide_diagnostics, only: point_diag

  implicit none

  private
  public :: glissade_calving_mask_init, glissade_calve_ice, &
            glissade_redistribute_unprotected_ice, &
            glissade_remove_icebergs, glissade_remove_isthmuses, &
            glissade_cull_calving_front, glissade_limit_cliffs,  &
            glissade_stress_tensor_eigenvalues, glissade_strain_rate_tensor_eigenvalues, &
            glissade_extrapolate_to_calving_front
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
       call glissade_get_masks(&
            nx,            ny,             &
            parallel,                      &
            thck,          topg,           &
            eus,           thklim,         &
            ice_mask,                      &
            ocean_mask = ocean_mask)

       ! Set the calving mask to include all ice-free ocean cells.
       ! Any ice entering these cells during the run will calve.

       ! Note: We tested an exception for cells with observed nonzero velocity
       ! (and hence ice present at the time of velocity observations) at vertices.
       ! The goal was to better represent the Thwaites shelf region, but the spin-up
       !  did not improve.
       ! Leaving the commented-out code in case we want to add something similar later.

       do j = 2, ny-1
          do i = 2, nx-1
             if (ocean_mask(i,j) == 1) then
!                if (usfc_obs(i-1,j)**2   + vsfc_obs(i-1,j)**2   > 0.0d0 .and. &
!                    usfc_obs(i,j)**2     + vsfc_obs(i,j)**2     > 0.0d0 .and. &
!                    usfc_obs(i-1,j-1)**2 + vsfc_obs(i-1,j-1)**2 > 0.0d0 .and. &
!                    usfc_obs(i,j-1)**2   + vsfc_obs(i,j-1)**2   > 0.0d0) then
!                   calving_mask(i,j) = 0
!                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
!                   if (verbose_calving) then   ! debug
!                      print*, 'ocean cell with uobs, vobs > 0: ig, jg =', iglobal, jglobal
!                   endif
!                else
!                   calving_mask(i,j) = 1   ! calve ice in this cell
!                endif
                calving_mask(i,j) = 1   ! calve ice in this cell
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
                                thck_pre_transport,      &  ! m
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
!    real(dp), intent(in)                     :: minthck             !> min thickness for ice at the calving front (m)
!    real(dp), intent(in)                     :: dthck_dx_cf         !> assumed max thickness gradient (m/m) at the subgrid CF
!    real(dp), dimension(:,:), intent(inout)  :: thck_effective      !> effective thickness for calving (m)
!    real(dp), dimension(:,:), intent(inout)  :: lateral_rate        !> lateral calving rate (m/s) at calving front
!    real(dp), dimension(:,:), intent(in)     :: tau_eigen1          !> first eigenvalue of 2D horizontal stress tensor (Pa)
!    real(dp), dimension(:,:), intent(in)     :: tau_eigen2          !> second eigenvalue of 2D horizontal stress tensor (Pa)
!    real(dp), dimension(:,:), intent(in)     :: eps_eigen1          !> first eigenvalue of 2D horizontal strain-rate tensor (1/s)
!    real(dp), dimension(:,:), intent(in)     :: eps_eigen2          !> second eigenvalue of 2D horizontal strain-rate tensor (1/s)
!    real(dp), dimension(:,:,:), intent(inout):: damage              !> 3D scalar damage parameter
!    real(dp), intent(in)                     :: damage_threshold    !> threshold value where ice is sufficiently damaged to calve
!    real(dp), intent(in)                     :: damage_constant     !> rate of change of damage (1/s) per unit stress (Pa)
!    integer,  dimension(:,:), intent(in)     :: calving_mask        !> integer mask: calve ice where calving_mask = 1
!    real(dp), dimension(:,:), intent(out)    :: calving_thck        !> thickness lost due to calving in each grid cell (m)

    integer, intent(in) :: itest, jtest, rtest                     !> coordinates of diagnostic point
    real(dp), intent(in)                      :: dt                !> model timestep (s)
    real(dp), intent(in)                      :: dx, dy            !> grid cell size in x and y directions (m)
    real(dp), dimension(:), intent(in)        :: sigma             !> vertical sigma coordinate
    real(dp), dimension(nx-1,ny-1), intent(in):: velnorm_mean      !> mean ice speed at vertices (m/s)
    real(dp), dimension(nx,ny), intent(in)    :: thck_pre_transport!> ice thickness (m) before doing transport, SMB, and BMB
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

    if (which_calving == CALVING_THCK_THRESHOLD) then

       call thickness_based_calving(&
            nx,                 ny,                    &
            dx,                 dy,                    &
            dt,                                        &
            itest,   jtest,     rtest,                 &
            parallel,                                  &
            which_ho_calving_front,                    &
            thck,               topg,                  &  ! m
            eus,                thklim,                &  ! m
            calving%thck_effective,                    &  ! m
            calving%minthck,                           &  ! m
            calving%dthck_dx_cf,                       &  ! m
            calving%timescale,                         &  ! s
            calving%calving_thck)                         ! m

    elseif (which_calving == EIGENCALVING) then

       ! Notes on eigencalving and damage-based calving:
       ! For both damage-based calving and eigencalving, the calving rate depends on eigenvalues of the 2D stress tensor.
       ! The main difference is that for eigencalving, the calving rate is based on current stresses
       !  at the calving front, whereas for damage-based calving, the calving rate is based on cumulative damage,
       !  which is generated in floating cells due to stresses and then is advected downstream.

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
            calving%thck_effective,                    &  ! m
!!            calving%minthck,                           &  ! m
            calving%dthck_dx_cf,                       &  ! m
            calving%lateral_rate,                      &  ! m/s
            calving%calving_thck)                         ! m

       ! Call thickness-based calving to remove thin floating ice

       if (calving%minthck > 0.0d0) then

          call thickness_based_calving(&
               nx,                 ny,                    &
               dx,                 dy,                    &
               dt,                                        &
               itest,   jtest,     rtest,                 &
               parallel,                                  &
               which_ho_calving_front,                    &
               thck,               topg,                  &  ! m
               eus,                thklim,                &  ! m
               calving%thck_effective,                    &  ! m
               calving%minthck,                           &  ! m
               calving%dthck_dx_cf,                       &  ! m
               calving%timescale,                         &  ! s
               calving%calving_thck)                         ! m

       endif

    elseif (which_calving == CALVING_DAMAGE) then

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
            calving%thck_effective,                    &  ! m
!!            calving%minthck,                           &  ! m
            calving%dthck_dx_cf,                       &  ! m
            calving%damage,                            &  !
            calving%lateral_rate,                      &  ! m/s
            calving%calving_thck)                         ! m

       ! Call thickness-based calving to remove thin floating ice

       if (calving%minthck > 0.0d0) then

          call thickness_based_calving(&
               nx,                 ny,                    &
               dx,                 dy,                    &
               dt,                                        &
               itest,   jtest,     rtest,                 &
               parallel,                                  &
               which_ho_calving_front,                    &
               thck,               topg,                  &  ! m
               eus,                thklim,                &  ! m
               calving%thck_effective,                    &  ! m
               calving%minthck,                           &  ! m
               calving%dthck_dx_cf,                       &  ! m
               calving%timescale,                         &  ! s
               calving%calving_thck)                         ! m

       endif

    elseif (which_calving == CF_ADVANCE_RETREAT_RATE) then

       call calving_front_advance_retreat(&
               nx,                 ny,                    &
               dx,                 dy,                    &
               dt,                                        &
               itest,   jtest,     rtest,                 &
               parallel,                                  &
               which_ho_calving_front,                    &
               thck_pre_transport,                        &  ! m
               thck,               topg,                  &  ! m
               eus,                thklim,                &  ! m
               calving%thck_effective,                    &  ! m
               calving%dthck_dx_cf,                       &  ! m
               calving%calving_thck)                         ! m

    else   ! other calving options
           !TODO - Put these in a separate subroutine

       ! Get masks.
       ! Use thickness limit of 0.0 instead of thklim so as to remove ice from any cell
       !  that meets the calving criteria, not just dynamically active ice.

       call glissade_get_masks(&
            nx,            ny,             &
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

  subroutine glissade_redistribute_unprotected_ice(&
       nx,              ny,       &
       itest,  jtest,   rtest,    &
       parallel,                  &
       protected_mask,            &
       thck_old,                  &
       thck_effective,            &
       thck,                      &
       calving_thck)

    use glide_diagnostics, only: point_diag

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    integer, dimension(nx,ny), intent(in) :: &
         protected_mask            ! mask to identify cells protected from ice removal or redistribution;
                                   ! includes land cells, full cells, and partial CF cells

    real(dp), dimension(nx,ny), intent(in) :: &
         thck_old,               & ! ice thickness (m) at start of timestep, before transport
         thck_effective            ! effective ice thickness (m) before transport;
                                   ! thck_effective < thck_pre_transport in partial CF cells

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m) before and after redistribution
         calving_thck              ! thickness (m) calved from each cell

    ! local variables

    integer :: i, j, ii, jj
    real(dp) :: thck_max           ! max thickness (m) in protected upstream cell
    real(dp) :: dthck              ! ice thickness (m) to be redistributed

    ! Write the initial thicknesses
    if (verbose_calving) then
       call point_diag(thck, 'Before redistribution, thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(thck_old, 'thck_old (m)', itest, jtest, rtest, 7, 7)
       call point_diag(thck_effective, 'thck_effective (m)', itest, jtest, rtest, 7, 7)
       call point_diag(protected_mask, 'protected_mask', itest, jtest, rtest, 7, 7)
    endif

    ! Identify unprotected ice with nonzero thickness.
    ! Instead of calving this ice, move it to a protected upstream cell
    ! (from which most or all of the ice likely arrived during transport).
    !TODO: Modify the logic to redistribute ice to more than one upstream cell.

    do j = 2, ny-1
       do i = 2, nx-1

          if (thck(i,j) > 0.0d0 .and. protected_mask(i,j) == 0) then

             ! Identify the protected edge neighbor with the thickest ice
             thck_max = 0
             if (thck_old(i-1,j) > thck_max .and. protected_mask(i-1,j) == 1) then
                thck_max = thck_old(i-1,j); ii = i-1; jj = j
             endif
             if (thck_old(i+1,j) > thck_max .and. protected_mask(i+1,j) == 1) then
                thck_max = thck_old(i+1,j); ii = i+1; jj = j
             endif
             if (thck_old(i,j-1) > thck_max .and. protected_mask(i,j-1) == 1) then
                thck_max = thck_old(i,j-1); ii = i; jj = j-1
             endif
             if (thck_old(i,j+1) > thck_max .and. protected_mask(i,j+1) == 1) then
                thck_max = thck_old(i,j+1); ii = i; jj = j+1
             endif

             ! Move ice from the unprotected cell to the edge neighbor
             ! (but not to exceed thck = thck_effective in the neighbor).
             ! If any ice remains in the unprotected cell, calve it now.
             if (thck_max > 0.0d0) then
                if (thck(ii,jj) + thck(i,j) < thck_effective(ii,jj)) then
                   dthck = thck(i,j)
                   thck(ii,jj) = thck(ii,jj) + dthck
                   thck(i,j) = 0.0d0
                else   ! fill the upstream cell up to H = H_eff, and calve the rest
                   dthck = thck_effective(ii,jj) - thck(ii,jj)
                   thck(ii,jj) = thck(ii,jj) + dthck
                   calving_thck(i,j) = calving_thck(i,j) + thck(i,j) - dthck
                   thck(i,j) = 0.0d0
                endif
             endif

          endif   ! unprotected and thck > 0
       enddo   ! i
    enddo   ! j

    call parallel_halo(thck, parallel)

    ! Write the final thicknesses
    if (verbose_calving) then
       call point_diag(thck, 'After redistribution, thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(calving_thck, 'calving_thck (m)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine glissade_redistribute_unprotected_ice

!---------------------------------------------------------------------------

  subroutine thickness_based_calving(&
       nx,                 ny,             &
       dx,                 dy,             &
       dt,                                 &
       itest,   jtest,     rtest,          &
       parallel,                           &
       which_ho_calving_front,             &
       thck,               topg,           &  ! m
       eus,                thklim,         &  ! m
       thck_effective,                     &  ! m
       calving_minthck,                    &  ! m
       dthck_dx_cf,                        &  ! m/m
       calving_timescale,                  &  ! s
       calving_thck)                          ! m

    ! Calve ice that is thinner than a prescribed threshold.
    ! This option requires a subgrid calving front scheme.
    ! The CF thinning rate is based not on the nominal ice thickness (H = thck),
    !  but on the effective calving thickness (thck_effective = H_eff),
    !  which depends on the max thickness of the edge-adjacent floating cells.

    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask
    use glide_diagnostics, only: point_diag

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
         topg                      ! bed topography (m)

    real(dp), intent(in) ::  &
         eus,                    & ! eustatic sea level (m)
         thklim                    ! minimum thickness for dynamically active grounded ice (m)

    real(dp), intent(in) :: &
         calving_minthck,        & ! min effective thickness (m) for CF cells
         dthck_dx_cf,            & ! assumed max thickness gradient (m/m) at the subgrid CF
         calving_timescale         ! timescale (s) for calving

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m)
         thck_effective,         & ! effective thickness for calving (m)
         calving_thck              ! thickness lost due to calving (m)

    ! local variables

    integer :: i, j

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
         thck_old,               & ! old value of thck (m)
         dt_calving,             & ! time remaining for calving (s)
         dt_calving_old            ! old value of dt_calving

    real(dp) :: &
         thinning_rate,          & ! vertical thinning rate (m/s)
         dthck                     ! thickness change (m)

    logical :: iterate_calving     ! Iterate the calving until dt_calving = 0 for all cells
    integer :: count

    ! Get masks

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
         dthck_dx_cf = dthck_dx_cf,         &
         dx = dx,       dy = dy,            &
         thck_effective = thck_effective,   &
         partial_cf_mask = partial_cf_mask, &
         full_mask = full_mask)

    if (verbose_calving) then
       call point_diag(calving_front_mask, 'Thickness-based calving, calving_front_mask', itest, jtest, rtest, 7, 7)
       call point_diag(thck, 'thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(thck_effective, 'thck_effective (m)', itest, jtest, rtest, 7, 7)
    endif

    ! Supplement the calving_front mask to include ocean cells with very thin ice (0 < thck > thklim).
    ! Thickness-based calving in these dynamically inactive cells reduces thickness oscillations around thklim.

    where (ocean_mask == 1 .and. thck_effective > 0.0d0)
       calving_front_mask = 1
    endwhere

    if (verbose_calving) then
       call point_diag(calving_front_mask, 'Adjusted calving_front_mask', itest, jtest, rtest, 7, 7)
    endif

    ! Setup for iteration.
    ! If a CF cell is removed entirely and time remains for calving,
    !  then recompute the CF and continue calving until dt_calving = 0.

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

       ! Apply thinning in calving-front cells whose effective thickness (H_eff = thck_effective)
       !  is less than a prescribed minimum value (Hc_min = calving_minthck).
       !
       ! The effective thinning rate is given by
       !
       !    dH_eff/dt = -(Hc_min - H_eff) / tau_c  where Hc_min > H_e
       !    dH_eff/dt = 0 elsewhere
       !
       ! where tau_c = calving%timescale.
       !
       ! The thinning rate applied to the mean cell thickness (thck) is given by
       !
       !    dH/dt = min(H/H_eff, 1) * dH_eff/dt
       !
       ! Thus, any ice with H_eff < Hc_min is removed on a time scale of tau_c.
       !TODO- Check the logic above. Is this what we do here? I don't see min(H/H_eff)

       do j = 2, ny-1
          do i = 2, nx-1
             if (calving_front_mask(i,j) == 1 .and. thck_effective(i,j) < calving_minthck) then

                if (calving_timescale > 0.0d0) then
                   thinning_rate = (calving_minthck - thck_effective(i,j)) / calving_timescale
                   dthck = thinning_rate * dt_calving(i,j)
                else
                   dthck = thck(i,j)
                endif

!                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
!                   write(6,*) ' '
!                   write(6,*) 'Thinning: r, i, j =', rtest, itest, jtest
!                   write(6,*) 'thck:', thck(i,j)
!                   write(6,*) 'thck_effective (m) =', thck_effective(i,j)
!                   write(6,*) 'calving_minthck (m) =', calving_minthck
!                   write(6,*) 'thinning rate (m/yr) =', thinning_rate * scyr
!                   write(6,*) 'dthck (m) =', dthck
!                endif

                if (dthck > thck(i,j)) then
                   ! calve the full column and compute the time remaining for more calving
                   dt_calving(i,j) = dt_calving(i,j) * (1.0d0 - thck(i,j)/dthck)
                   iterate_calving = .true.
                   if (verbose_calving .and. this_rank == rtest) then
                      write(6,*) 'calving time remaining (yr):', i, j, dt_calving(i,j)/scyr
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

             endif   ! calving_front cell with H_eff < Hc_min

          enddo   ! i
       enddo   ! j

       call parallel_halo(thck, parallel)
       call parallel_halo(calving_thck, parallel)

       ! Check whether any processors have iterate_calving = T.
       ! If so, then continue.  If not, then exit.
       iterate_calving = parallel_reduce_log_or(iterate_calving)

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

          call glissade_calving_front_mask(&
               nx,            ny,                 &
               which_ho_calving_front,            &
               parallel,                          &
               thck,          topg,               &
               eus,                               &
               ice_mask,      floating_mask,      &
               ocean_mask,    land_mask,          &
               calving_front_mask)

          where (ocean_mask == 1 .and. thck_effective > 0.0d0)
             calving_front_mask = 1
          endwhere

          ! Save the old dt_calving values and compute new values
          ! Loop through cells; if a new CF cell has ice-free neighbors with dt_calving > 0,
          !  then use this time to calve.  If it has multiple such neighbors, then use the largest value.

          dt_calving_old = dt_calving
          dt_calving = 0.0d0

          do j = 2, ny-1
             do i = 2, nx-1
                if (calving_front_mask_old(i,j) == 0 .and. calving_front_mask(i,j) == 1) then
                   dt_calving(i,j) = max(dt_calving_old(i-1,j), dt_calving_old(i+1,j), &
                                         dt_calving_old(i,j-1), dt_calving_old(i,j+1))
                   if (verbose_calving .and. this_rank == rtest .and. dt_calving(i,j) > 0.0d0) then
                      write(6,*) 'rank, i, j, new dt_calving (yr):', this_rank, i, j, dt_calving(i,j)/scyr
                   endif
                endif
             enddo
          enddo

       endif   ! iterate_calving

    enddo   ! do while iterate_calving = T

  end subroutine thickness_based_calving

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
       thck_effective,                    &
       dthck_dx_cf,                       &
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
         dthck_dx_cf               ! assumed max thickness gradient (m/m) at the subgrid CF

    real(dp), dimension(nx,ny), intent(out) :: &
         lateral_rate              ! lateral rate of calving (m/s)

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m)
         thck_effective,         & ! effective thickness for calving (m)
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
         dthck_dx_cf = dthck_dx_cf,         &
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

       if (verbose_calving) then
          call point_diag(thck_effective, 'thck_effective (m)', itest, jtest, rtest, 7, 7)
          call point_diag(thck, 'thck (m)', itest, jtest, rtest, 7, 7)
          call point_diag(partial_cf_mask, 'partial_cf_mask', itest, jtest, rtest, 7, 7)
          call point_diag(full_mask, 'full_mask', itest, jtest, rtest, 7, 7)
          call point_diag(calving_front_mask, 'calving_front_mask', itest, jtest, rtest, 7, 7)
       endif

       ! Main calving loop
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (calving_front_mask(i,j) == 1 .and. dt_calving(i,j) > 0.0d0) then

                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank == rtest .and. count > 1) then
                   write(6,*) 'Next round of calving: count, i, j, dt_calving (yr):', &
                        count, i, j, dt_calving(i,j)/scyr
                endif

                ! Compute a rate factor to account for the length of the calving front.
                ! A typical CF cell will have ice-free ocean on just one edge.
                ! Cells with 2 or 3 ocean neighbors have a longer CF and will calve more quickly.

                ocean_neighbor_count = &
                     ocean_mask(i-1,j) + ocean_mask(i+1,j) + ocean_mask(i,j-1) + ocean_mask(i,j+1)

                lateral_rate_factor = ocean_neighbor_count

                if (verbose_calving .and. this_rank == rtest .and. abs(i-itest) <= 3 .and. abs(j-jtest) <= 3) then
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

       ! Check whether any processors have iterate_calving = T.
       ! If so, then continue.  If not, then exit.
       iterate_calving = parallel_reduce_log_or(iterate_calving)

       if (iterate_calving) then

          ! Recompute the CF mask
          ! This could result in new CF cells that were previously in the interior.
          ! Note: thck_effective keeps the values computed above.

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

          call glissade_calving_front_mask(&
               nx,            ny,                 &
               which_ho_calving_front,            &
               parallel,                          &
               thck,          topg,               &
               eus,                               &
               ice_mask,      floating_mask,      &
               ocean_mask,    land_mask,          &
               calving_front_mask)

          ! Save the old dt_calving values and compute new values
          ! Loop through cells; if a new CF cell has ice-free neighbors with dt_calving > 0,
          !  then use this time to calve.  If it has multiple such neighbors, then use the largest value.

          dt_calving_old = dt_calving
          dt_calving = 0.0d0

          do j = 2, ny-1
             do i = 2, nx-1
                if (calving_front_mask_old(i,j) == 0 .and. calving_front_mask(i,j) == 1) then
                   dt_calving(i,j) = max(dt_calving_old(i-1,j), dt_calving_old(i+1,j), &
                                         dt_calving_old(i,j-1), dt_calving_old(i,j+1))
                   iterate_calving = .true.
                   if (verbose_calving .and. i==itest .and. j==jtest .and. &
                        this_rank == rtest .and. dt_calving(i,j) > 0.0d0) then
                      write(6,*) 'rank, i, j, new dt_calving (yr):', this_rank, i, j, dt_calving(i,j)/scyr
                   endif
                endif
             enddo
          enddo

       endif   ! iterate_calving

    enddo   ! do while iterate_calving = T

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
       damage_threshold,                          &
       thck_effective,                            &
       dthck_dx_cf,                               &
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
         dthck_dx_cf               ! assumed max thickness gradient (m/m) at the subgrid CF

    real(dp), dimension(nz-1,nx,ny), intent(inout) :: &
         damage                    ! 3D damage tracer, 0 > damage < 1

    real(dp), dimension(nx,ny), intent(out) :: &
         lateral_rate              ! lateral rate of calving (m/s)

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m)
         thck_effective,         & ! effective thickness for calving (m)
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
!    call parallel_halo(eps_eigen1, parallel)
!    call parallel_halo(eps_eigen2, parallel)

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
         dthck_dx_cf = dthck_dx_cf,         &
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
                ! (1) d_damage_dt is proportional to the weighted stress eigenvalues
                ! (2) d_damage_dt is proportional to the weighted stress eigenvalues, divided by H
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
          call point_diag(thck_effective, 'thck_effective (m)', itest, jtest, rtest, 7, 7)
          call point_diag(thck, 'thck (m)', itest, jtest, rtest, 7, 7)
          call point_diag(partial_cf_mask, 'partial_cf_mask', itest, jtest, rtest, 7, 7)
          call point_diag(full_mask, 'full_mask', itest, jtest, rtest, 7, 7)
          call point_diag(calving_front_mask, 'calving_front_mask', itest, jtest, rtest, 7, 7)
       endif

       ! Main calving loop
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (calving_front_mask(i,j) == 1 .and. damage_column(i,j) > damage_threshold  &
                  .and. dt_calving(i,j) > 0.0d0) then

                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank == rtest .and. count > 1) then
                   write(6,*) 'Next round of calving: count, rank, i, j, dt_calving (yr):', &
                        count, this_rank, i, j, dt_calving(i,j)/scyr
                endif

                ! Compute a damage fraction, increasing linearly from 0 (at a damage threshold) to 1.
                damage_frac = (damage_column(i,j) - damage_threshold) / (1.0d0 - damage_threshold)

                ! Compute a rate factor to account for the length of the calving front.
                ! A typical CF cell will have ice-free ocean on just one edge.
                ! Cells with 2 or 3 ocean neighbors have a longer CF and will calve more quickly.

                ocean_neighbor_count = &
                     ocean_mask(i-1,j) + ocean_mask(i+1,j) + ocean_mask(i,j-1) + ocean_mask(i,j+1)

                lateral_rate_factor = ocean_neighbor_count

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

       ! Check whether any processors have iterate_calving = T.
       ! If so, then continue.  If not, then exit.
       iterate_calving = parallel_reduce_log_or(iterate_calving)

       if (iterate_calving) then

          ! Recompute the CF mask
          ! This could result in new CF cells that were previously in the interior.
          ! Note: thck_effective keeps the values computed above.

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

          call glissade_calving_front_mask(&
               nx,            ny,                 &
               which_ho_calving_front,            &
               parallel,                          &
               thck,          topg,               &
               eus,                               &
               ice_mask,      floating_mask,      &
               ocean_mask,    land_mask,          &
               calving_front_mask,                &
               dthck_dx_cf = dthck_dx_cf,         &
               dx = dx,       dy = dy,            &
               thck_effective = thck_effective,   &
               partial_cf_mask = partial_cf_mask, &
               full_mask = full_mask)

          ! Save the old dt_calving values and compute new values
          ! Loop through cells; if a new CF cell has ice-free neighbors with dt_calving > 0,
          !  then use this time to calve.  If it has multiple such neighbors, then use the largest value.

          dt_calving_old = dt_calving
          dt_calving = 0.0d0

          do j = 2, ny-1
             do i = 2, nx-1
                if (calving_front_mask_old(i,j) == 0 .and. calving_front_mask(i,j) == 1  &
                     .and. damage_column(i,j) > damage_threshold) then   ! new CF cell that can calve
                   dt_calving(i,j) = max(dt_calving_old(i-1,j), dt_calving_old(i+1,j), &
                                         dt_calving_old(i,j-1), dt_calving_old(i,j+1))
                   if (verbose_calving .and. this_rank == rtest .and. dt_calving(i,j) > 0.0d0) then
                      write(6,*) 'rank, i, j, new dt_calving (yr):', this_rank, i, j, dt_calving(i,j)/scyr
                   endif
                endif
             enddo
          enddo

       endif   ! iterate_calving

    enddo   ! do while iterate_calving = T

  end subroutine damage_based_calving

!---------------------------------------------------------------------------

  subroutine calving_front_advance_retreat(&
       nx,                 ny,             &
       dx,                 dy,             &
       dt,                                 &
       itest,   jtest,     rtest,          &
       parallel,                           &
       which_ho_calving_front,             &
       thck_pre_transport,                 &  ! m
       thck,               topg,           &  ! m
       eus,                thklim,         &  ! m
       thck_effective,                     &  ! m
       dthck_dx_cf,                        &  ! m/m
       calving_thck)                          ! m

    ! Calve ice based on a prescribed rate of calving front advance or retreat.
    ! This option requires a subgrid calving front scheme.

    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask
    use glide_diagnostics, only: point_diag

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
         thck_pre_transport,     & ! ice thickness (m) before doing transport, SMB, and BMB;
                                   ! determines the calving needed to keep the CF stationary
         topg                      ! bed topography (m)

    real(dp), intent(in) ::  &
         eus,                    & ! eustatic sea level (m)
         thklim                    ! minimum thickness for dynamically active grounded ice (m)

    real(dp), intent(in) :: &
         dthck_dx_cf               ! assumed max thickness gradient (m/m) at the subgrid CF

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m)
         thck_effective,         & ! effective thickness for calving (m)
         calving_thck              ! thickness lost due to calving (m)

    ! local variables

    integer :: i, j

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
         thck_old,               & ! old value of thck (m)
         dt_calving_old,         & ! old value of dt_calving
         dt_calving                ! time remaining for calving (s)

    real(dp) :: &
         thinning_rate,          & ! vertical thinning rate (m/s)
         dthck                     ! thickness change (m)

    !WHL - One might ask, why iterate? Here is one situation in which it matters.
    !      Consider a setting like MISMIP+, with flow mainly in the x-direction.
    !      We start with a sharp calving front, with thick ice to the left and ice-free ocean to the right.
    !      During transport (which precedes calving), ice flows into the first row of ice-free cells.
    !      Suppose we are prescribing CF retreat. This means that all the ice in the partial cells
    !       should be calved, and in addition, the full cells that were in the interior should be thinned.
    !      If we calve ice in the partial cells and do nothing more, the CF will not retreat at the desired speed.

    logical :: iterate_calving     ! Iterate the calving until dt_calving = 0 for all cells
    integer :: count

    !WHL - temporary
    real(dp) :: cf_advance_retreat_rate  ! prescribed rate of calving-front advance or retreat (m/s);
                                         ! positive for advance

    real(dp), dimension(nx,ny) :: &
         dthck_dt_transport              ! rate of thickness change from transport, SMB and BMB

    ! Get masks

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
         dthck_dx_cf = dthck_dx_cf,         &
         dx = dx,       dy = dy,            &
         thck_effective = thck_effective,   &
         partial_cf_mask = partial_cf_mask, &
         full_mask = full_mask)

    !WHL - For now, compute the CF advance/retreat rate locally.
    !      Later, pass it in

    cf_advance_retreat_rate =    0.0d0       ! stationary
!    cf_advance_retreat_rate = -300.d0/scyr   ! retreat of 300 m/yr, converted to m/s
!    cf_advance_retreat_rate =  300.d0/scyr   ! advance of 300 m/yr, converted to m/s

    dthck_dt_transport = (thck - thck_pre_transport) / dt

    if (verbose_calving) then
       call point_diag(thck, 'Prescribe CF advance/retreat: thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(thck_pre_transport, 'thck_pre_transport (m)', itest, jtest, rtest, 7, 7)
       call point_diag(dthck_dt_transport*scyr, 'dH/dt transport (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(thck_effective, 'thck_effective (m)', itest, jtest, rtest, 7, 7)
       call point_diag(calving_front_mask, 'calving_front_mask', itest, jtest, rtest, 7, 7)
       call point_diag(partial_cf_mask, 'partial_cf_mask', itest, jtest, rtest, 7, 7)
       call point_diag(full_mask, 'full_mask', itest, jtest, rtest, 7, 7)
    endif

    ! Loop over locally owned cells
    ! Calving occurs only in CF cells: cells with one or more ocean neighbors.
    ! If time remains after one or more columns have calved completely, iterate until dt_calving = 0 everywhere.

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
          call point_diag(dt_calving/scyr, 'dt_calving (yr)', itest, jtest, rtest, 7, 7)
       endif

       iterate_calving = .false.
       thck_old = thck   ! save current value of thck

       !TODO - For iterations 2 and higher, simply apply the remaining dthck to the upstream cell?
       ! Main calving loop
       do j = 2, ny-1
          do i = 2, nx-1

             if (calving_front_mask(i,j) == 1 .and. dt_calving(i,j) > 0.0d0) then

                ! First, compute the calving needed to compensate for the thickness change in partial CF cells.
                ! This step should hold the calving front stationary by balancing transport/SMB/BMB with calving.

                dthck = max(dthck_dt_transport(i,j), 0.0d0) * dt_calving(i,j)

                ! If the CF is supposed to retreat, then increase the calving.
                ! If the CF is supposed to advance, then reduce the calving.
                ! Note: Some calving might already have been done in an unprotected cell just past the current CF.
                !       If so, it is allowed to undo this calving so the CF can advance.

                if (cf_advance_retreat_rate /= 0.0d0) then

                   ! increase dthck if the CF is retreating, decrease if the CF is advancing
                   dthck = dthck - (cf_advance_retreat_rate * dt_calving(i,j)) * thck_effective(i,j) / sqrt(dx*dy)

                   ! Make sure the sum of calving_thck (as previously computed) and dthck is non-negative
                   ! This prevents unphysical negative calving and limits the rate of CF advance.
                   if (calving_thck(i,j) + dthck < 0.0d0) dthck = -calving_thck(i,j)

                endif


                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                   write(6,*) ' '
                   write(6,*) 'count, rank, i, j, dt_calving (yr):', &
                        count, this_rank, i, j, dt_calving(i,j)/scyr
                   write(6,*) 'cf_advance_retreat_rate (m/yr) =', cf_advance_retreat_rate * scyr
                   write(6,*) 'thck:', thck(i,j)
                   write(6,*) 'thck_effective (m) =', thck_effective(i,j)
                   write(6,*) 'dthck transport (m) =', dthck_dt_transport(i,j) * dt_calving(i,j)
                   write(6,*) 'dthck calving (m)     =', dthck
                endif

                ! Apply this thinning to CF cells

                if (dthck > thck(i,j)) then
                   ! calve the full column and compute the time remaining for more calving
                   dt_calving(i,j) = dt_calving(i,j) * (1.0d0 - thck(i,j)/dthck)
                   iterate_calving = .true.
                   if (verbose_calving .and. this_rank == rtest) then
                      write(6,*) 'thck, dthck, ratio:', thck(i,j), dthck, thck(i,j)/dthck
                      write(6,*) 'calving time remaining (yr):', i, j, dt_calving(i,j)/scyr
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

             endif   ! calving_front cell

          enddo   ! i
       enddo   ! j

       call parallel_halo(thck, parallel)
       call parallel_halo(calving_thck, parallel)

       ! Check whether any processors have iterate_calving = T.
       ! If so, then continue.  If not, then exit.
       iterate_calving = parallel_reduce_log_or(iterate_calving)

       if (iterate_calving) then

          ! Recompute the CF mask
          ! This could result in new CF cells that were previously in the interior.
          ! Note: thck_effective keeps the values computed above.

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

          call glissade_calving_front_mask(&
               nx,            ny,                 &
               which_ho_calving_front,            &
               parallel,                          &
               thck,          topg,               &
               eus,                               &
               ice_mask,      floating_mask,      &
               ocean_mask,    land_mask,          &
               calving_front_mask,                &
               dthck_dx_cf = dthck_dx_cf,         &
               dx = dx,       dy = dy,            &
               thck_effective = thck_effective,   &
               partial_cf_mask = partial_cf_mask, &
               full_mask = full_mask)

          ! Save the old dt_calving values and compute new values
          ! Loop through cells; if a new CF cell has ice-free neighbors with dt_calving > 0,
          !  then use this time to calve.  If it has multiple such neighbors, then use the largest value.

          dt_calving_old = dt_calving
          dt_calving = 0.0d0

          do j = 2, ny-1
             do i = 2, nx-1
                if (calving_front_mask_old(i,j) == 0 .and. calving_front_mask(i,j) == 1) then   ! new CF cell that can calve
                   dt_calving(i,j) = max(dt_calving_old(i-1,j), dt_calving_old(i+1,j), &
                                         dt_calving_old(i,j-1), dt_calving_old(i,j+1))
                   if (verbose_calving .and. i==itest .and. j==jtest .and. &
                        this_rank == rtest .and. dt_calving(i,j) > 0.0d0) then
                      write(6,*) 'rank, i, j, new dt_calving (yr):', this_rank, i, j, dt_calving(i,j)/scyr
                   endif
                endif
             enddo
          enddo

       endif   ! iterate_calving

    enddo   ! do while iterate_calving = T

  end subroutine calving_front_advance_retreat

!---------------------------------------------------------------------------

  !TODO - Remove this subroutine? Ideally, unstable ice can be removed by other means.
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

    do n = 1, ncull_calving_front

       ! calculate masks
       ! Note: Passing in thklim = 0.0 does not work because it erroneously counts thin floating cells as active.
       !       Then the algorithm can fail to identify floating regions that should be removed
       !       (since they are separated from any active cells).

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
            calving_front_mask)

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
       f_ground_threshold,          &
       thck,                        &
       f_ground_cell,               &
       ice_mask,                    &
       land_mask,                   &
       calving_thck)

    ! Remove any icebergs. 
        
    ! The algorithm is as follows:
    ! (1) Mark all cells with ice  with the initial color.
    !     Mark other cells with the boundary color.
    ! (2) Seed the fill by giving grounded ice cells the fill color.
    ! (3) Recursively fill all cells that are connected to filled cells by a path
    !     that passes through ice-covered cells only.
    ! (4) Repeat the recursion as necessary to spread the fill to adjacent processors.
    ! (5) Once the fill is done, any floating cells that still have the initial color
    !     are considered to be icebergs and are removed.
    !
    ! Notes:
    ! (1) Grounded cells must have f_ground_cell > f_ground_threshold to seed the fill.
    ! (2) The recursive fill applies to edge neighbors, not corner neighbors.
    !     The path back to grounded ice must go through edges, not corners.
    ! (3) Should have thklim > 0.  With a limit of 0.0, very thin floating cells
    !     can be wrongly counted as active, and icebergs can be missed.
    ! (4) Land-based cells that still have the initial color are not marked as icebergs.

    use glissade_masks, only: glissade_fill_with_buffer, initial_color, fill_color, boundary_color

    integer, intent(in) :: nx, ny                                !> horizontal grid dimensions
    type(parallel_type), intent(in) :: parallel                  !> info for parallel communication
    integer, intent(in) :: itest, jtest, rtest                   !> coordinates of diagnostic point
    real(dp), intent(in) :: f_ground_threshold                   !> threshold for counting cells as grounded

    real(dp), dimension(nx,ny), intent(inout) :: thck            !> ice thickness
    real(dp), dimension(nx,ny), intent(in)    :: f_ground_cell   !> grounded fraction in each grid cell
    integer,  dimension(nx,ny), intent(in)    :: ice_mask        !> = 1 where ice is present (thck > thklim), else = 0
    integer,  dimension(nx,ny), intent(in)    :: land_mask       !> = 1 where topg - eus >= 0, else = 0
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

    if (verbose_calving) then
       call point_diag(thck, 'Remove icebergs, thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(ice_mask, 'ice_mask', itest, jtest, rtest, 7, 7)
       call point_diag(f_ground_cell, 'f_ground_cell', itest, jtest, rtest, 7, 7)
    endif
    
    ! Initialize iceberg removal
    ! Note: Any cell with ice receives the initial color.
    !       TODO - The comments below no longer apply?
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

    ! Loop through cells, identifying cells with grounded ice.
    ! Fill each grounded cell and then recursively fill neighbor cells, whether grounded or not.
    ! We may have to do this several times to incorporate connections between neighboring processors.

    max_iter = max(parallel%ewtasks, parallel%nstasks)
    global_count_save = 0

    do iter = 1, max_iter

       if (iter == 1) then   ! identify grounded cells that can seed the fill

          do j = 1, ny
             do i = 1, nx

                ! Cells that are firmly grounded can seed the fill.
                ! It is somewhat arbitrary how to define "firmly grounded", but here we require
                !  that f_ground_cell exceeds a threshold value defined above.
                ! Note: If running without a GLP, then f_ground_cell is binary, either 0 or 1.

                if (ice_mask(i,j) == 1 .and. f_ground_cell(i,j) >= f_ground_threshold) then  ! grounded ice

                   if (color(i,j) /= boundary_color .and. color(i,j) /= fill_color) then

                      ! assign the fill color to this cell, and recursively fill neighbor cells
                      !TODO - Use glissade_fill instead of glissade_fill_with_buffer?  (Here and below)
                      call glissade_fill_with_buffer(&
                           nx,    ny,    &
                           i,     j,     &
                           color, ice_mask)

                   endif

                endif
             enddo
          enddo

       else  ! count > 1

          ! Check for halo cells that were just filled on neighbor processors
          ! TODO - Is the 'ice_mask = 1' logic redundant? (Here and below?)  Use glissade_fill without buffer?
          call parallel_halo(color, parallel)

          ! west halo layer
          i = nhalo
          do j = 1, ny
             if (color(i,j) == fill_color .and. ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(&
                     nx,    ny,    &
                     i+1,   j,     &
                     color, ice_mask)
             endif
          enddo

          ! east halo layers
          i = nx - nhalo + 1
          do j = 1, ny
             if (color(i,j) == fill_color .and. ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(&
                     nx,    ny,    &
                     i-1,   j,     &
                     color, ice_mask)
             endif
          enddo

          ! south halo layer
          j = nhalo
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(&
                     nx,    ny,    &
                     i,     j+1,   &
                     color, ice_mask)
             endif
          enddo

          ! north halo layer
          j = ny-nhalo+1
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color .and. ice_mask(i,j) == 1) then
                call glissade_fill_with_buffer(&
                     nx,    ny,    &
                     i,     j-1,   &
                     color, ice_mask)
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

    !TODO - Remove the exception?
    ! Icebergs are cells that still have the initial color and are not on land.
    ! Remove ice in these cells, adding it to the calving field.
    ! Note: There is an exception for cells that are
    !       (1) adjacent to at least one ice-covered cell (sharing an edge), and
!!    !       (2) connected diagonally to an active cell with the fill color.
    !       (2) connected diagonally to an ice-covered cell with the fill color.
    !       Such cells are considered part of the inactive calving front and are
    !        allowed to continue filling instead of calving.
    ! Allow land-based cells to be removed if f_ground < f_ground_threshold

    do j = 2, ny-1
       do i = 2, nx-1

          if (color(i,j) == initial_color .and. land_mask(i,j) == 0) then
             if (  ( color(i-1,j+1)==fill_color .and. ice_mask(i-1,j+1)==1 .and. &
                       (ice_mask(i-1,j)==1 .or. ice_mask(i,j+1)==1) ) &
              .or. ( color(i+1,j+1)==fill_color .and. ice_mask(i+1,j+1)==1 .and. &
                       (ice_mask(i+1,j)==1 .or. ice_mask(i,j+1)==1) ) &
              .or. ( color(i-1,j-1)==fill_color .and. ice_mask(i-1,j-1)==1 .and. &
                       (ice_mask(i-1,j)==1 .or. ice_mask(i,j-1)==1) ) &
              .or. ( color(i+1,j-1)==fill_color .and. ice_mask(i+1,j-1)==1 .and. &
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
       f_ground_threshold,          &
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
    real(dp), intent(in) :: f_ground_threshold          !> threshold for counting cells as grounded

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
    !  f_ground_threshold is used to identify weakly grounded cells.

    ! An isthmus cell has ice-free ocean or thin floating ice on each side:
    !  isthmus_thck_threshold is used to identify thin floating ice.
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
          if (floating_mask(i,j) == 1 .or. f_ground_cell(i,j) < f_ground_threshold) then
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
       taumax_cliff,   cliff_timescale, &
       thck,           topg,            &
       eus,            thklim,          &
       calving_thck)

    ! Impose a thickness limit on marine ice cliffs.
    ! These are defined as grounded marine-based cells adjacent to ice-free ocean.
    ! Ice removed from cliffs is added to the calving flux.

    use glissade_masks

    integer, intent(in)  :: nx, ny                      !> horizontal grid dimensions
    type(parallel_type), intent(in) :: parallel         !> info for parallel communication
    integer, intent(in)  :: itest, jtest, rtest         !> coordinates of diagnostic point
    real(dp), intent(in) :: dt                          !> model timestep (s)
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
         marine_cliff_mask     ! = 1 where ice is grounded and marine-based and borders at least one ocean cell

    real(dp), dimension(nx,ny) ::  &
         thckmax_cliff         ! max stable ice thickness in marine_cliff cells

    real(dp) :: &
         thinning_rate,        & ! vertical thinning rate (m/s)
         dthck,                & ! thickness change (m)
         factor

    ! Update masks, including the marine_cliff mask.

    call glissade_get_masks(&
         nx,            ny,             &
         parallel,                      &
         thck,          topg,           &
         eus,           thklim,         &
         ice_mask,                      &
         floating_mask = floating_mask, &
         ocean_mask = ocean_mask,       &
         land_mask = land_mask)

    call glissade_marine_cliff_mask(&
         nx,            ny,                &
         ice_mask,      floating_mask,     &
         land_mask,     ocean_mask,        &
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

  subroutine glissade_stress_tensor_eigenvalues(&
       nx,    ny,   nz,   &
       sigma,             &
       tau,               &
       tau_eigen1,        &
       tau_eigen2)

    ! Compute the eigenvalues of the 2D horizontal stress tensor.
    ! These are used for eigencalving and damage-based calving.

    use glimmer_paramets, only: tau0

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny, nz                ! grid dimensions

    real(dp), dimension(nz), intent(in) :: &
         sigma                     ! vertical sigma coordinate

    type(glide_tensor), intent(in) :: &
         tau                       ! 3D stress tensor (Pa)

    real(dp), dimension(nx,ny), intent(out) :: &
         tau_eigen1, tau_eigen2    ! eigenvalues of 2D horizontal stress tensor (Pa)

    ! local variables

    integer :: i, j, k
    real(dp) :: a, b, c, dsigma, root, lambda1, lambda2
    real(dp) :: tau_xx, tau_yy, tau_xy   ! vertically averaged stress tensor components

    tau_eigen1 = 0.0d0
    tau_eigen2 = 0.0d0

    do j = 1, ny
       do i = 1, nx

          ! compute vertically averaged stress components
          tau_xx = 0.0d0
          tau_yy = 0.0d0
          tau_xy = 0.0d0

          do k = 1, nz-1
             dsigma = sigma(k+1) - sigma(k)
             tau_xx = tau_xx + tau0 * tau%xx(k,i,j) * dsigma
             tau_yy = tau_yy + tau0 * tau%yy(k,i,j) * dsigma
             tau_xy = tau_xy + tau0 * tau%xy(k,i,j) * dsigma
          enddo

          ! compute the eigenvalues of the vertically integrated stress tensor
          a = 1.0d0
          b = -(tau_xx + tau_yy)
          c = tau_xx*tau_yy - tau_xy*tau_xy
          if (b*b - 4.0d0*a*c > 0.0d0) then   ! two real eigenvalues
             root = sqrt(b*b - 4.0d0*a*c)
             lambda1 = (-b + root) / (2.0d0*a)
             lambda2 = (-b - root) / (2.0d0*a)
             if (lambda1 > lambda2) then
                tau_eigen1(i,j) = lambda1
                tau_eigen2(i,j) = lambda2
             else
                tau_eigen1(i,j) = lambda2
                tau_eigen2(i,j) = lambda1
             endif
          endif  ! b^2 - 4ac > 0

       enddo   ! i
    enddo   ! j

  end subroutine glissade_stress_tensor_eigenvalues

!---------------------------------------------------------------------------

  subroutine glissade_strain_rate_tensor_eigenvalues(&
       nx,    ny,   nz,          &
       sigma,                    &
       strain_rate,              &
       eps_eigen1,  eps_eigen2,  &
       tau,         efvs,  &
       divu,        shear)

    ! Compute the eigenvalues of the 2D horizontal strain rate tensor.
    ! These can be used for eigencalving and damage-based calving, or for diagnostics.
    ! There are two ways to call the subroutine:
    ! (1) Pass in the strain rate tensor and compute the eigenvalues directly.
    ! (2) Pass in the stress tensor as an optional argument, compute the strain rate tensor
    !     from the stress tensor and effective viscosity, and then compute the eigenvalues.

    use glimmer_paramets, only: evs0, tau0

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny, nz                ! grid dimensions

    real(dp), dimension(nz), intent(in) :: &
         sigma                     ! vertical sigma coordinate

    type(glide_tensor), intent(inout) :: &
         strain_rate               ! 3D strain rate tensor
                                   ! intent(out) if computed from tau and efvs

    real(dp), dimension(nx,ny), intent(out) :: &
         eps_eigen1, eps_eigen2    ! eigenvalues of 2D horizontal stress tensor (1/s)

    type(glide_tensor), intent(in), optional :: &
         tau                       ! 3D stress tensor (Pa)

    real(dp), dimension(nz-1,nx,ny), intent(in), optional :: &
         efvs                      ! effective viscosity (Pa s)

    real(dp), dimension(nx,ny), intent(out), optional :: &
         divu,                   & ! divergence of horizontal flow (1/s)
         shear                     ! shear-related invariant of horizontal flow (1/s)
                                   ! not strictly shear since it includes a tensile term
    ! local variables

    integer :: i, j, k
    real(dp) :: a, b, c, dsigma, root, lambda1, lambda2
    real(dp) :: eps_xx, eps_yy, eps_xy   ! vertically averaged strain rate tensor components

    ! Optionally, compute the strain rate tensor from the stress tensor and effective viscosity

    if (present(tau) .and. present(efvs)) then

       where (efvs > 0.0d0)
          strain_rate%scalar = tau0 * tau%scalar / (2.d0 * evs0 * efvs)
          strain_rate%xz = tau0 * tau%xz / (2.d0 * evs0 * efvs)
          strain_rate%yz = tau0 * tau%yz / (2.d0 * evs0 * efvs)
          strain_rate%xx = tau0 * tau%xx / (2.d0 * evs0 * efvs)
          strain_rate%yy = tau0 * tau%yy / (2.d0 * evs0 * efvs)
          strain_rate%xy = tau0 * tau%xy / (2.d0 * evs0 * efvs)
       elsewhere
          strain_rate%scalar = 0.0d0
          strain_rate%xz = 0.0d0
          strain_rate%yz = 0.0d0
          strain_rate%xx = 0.0d0
          strain_rate%yy = 0.0d0
          strain_rate%xy = 0.0d0
       endwhere
    endif

    ! Compute the eigenvalues of the 2D horizontal strain rate tensor

    eps_eigen1 = 0.0d0
    eps_eigen2 = 0.0d0

    do j = 1, ny
       do i = 1, nx

          ! compute vertically averaged strain rate components
          eps_xx = 0.0d0
          eps_yy = 0.0d0
          eps_xy = 0.0d0

          do k = 1, nz-1
             dsigma = sigma(k+1) - sigma(k)
             eps_xx = eps_xx + strain_rate%xx(k,i,j) * dsigma
             eps_yy = eps_yy + strain_rate%yy(k,i,j) * dsigma
             eps_xy = eps_xy + strain_rate%xy(k,i,j) * dsigma
          enddo

          ! compute the eigenvalues of the vertically integrated strain rate tensor
          a = 1.0d0
          b = -(eps_xx + eps_yy)
          c = eps_xx*eps_yy - eps_xy*eps_xy
          if (b*b - 4.0d0*a*c > 0.0d0) then   ! two real eigenvalues
             root = sqrt(b*b - 4.0d0*a*c)
             lambda1 = (-b + root) / (2.0d0*a)
             lambda2 = (-b - root) / (2.0d0*a)
             if (lambda1 > lambda2) then
                eps_eigen1(i,j) = lambda1
                eps_eigen2(i,j) = lambda2
             else
                eps_eigen1(i,j) = lambda2
                eps_eigen2(i,j) = lambda1
             endif
          endif  ! b^2 - 4ac > 0

          ! Optionally, compute two other invariants of the horizontal flow:
          !    divu = eps_xx + eps_yy
          !    shear = sqrt{[(eps_xx - eps_yy)/2]^2 + eps_xy^2}
          ! These are related to the eigenvalues as:
          !    eps1 = divu + shear
          !    eps2 = divu - shear
          if (present(divu)) divu(i,j)  = (eps_xx + eps_yy)/2.0d0
          if (present(shear)) &
               shear(i,j) = sqrt(((eps_xx - eps_yy)/2.0d0)**2 + eps_xy**2)

       enddo   ! i
    enddo   ! j

  end subroutine glissade_strain_rate_tensor_eigenvalues

!---------------------------------------------------------------------------

  subroutine glissade_extrapolate_to_calving_front(&
       nx,              ny,            &
       partial_cf_mask, full_mask,     &
       field1,          field2)

    ! Extrapolate values from full cells to partial calving-front cells.
    ! This can be useful for dynamic fields whose values may be unrealistic
    !  in partly filled CF cells.

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny                ! horizontal grid dimensions

    integer, dimension(nx,ny), intent(in) :: &
         partial_cf_mask,    & ! = 1 for partial ice-covered calving-front cells, else = 0
         full_mask             ! = 1 for full ice-covered cells, else = 0

    real(dp), dimension(nx,ny), intent(out) :: &
         field1                ! field with values to be extrapolated to partial CF cells

    real(dp), dimension(nx,ny), intent(out), optional :: &
         field2                ! second field with values to be extrapolated to partial CF cells

    ! local variables

    integer :: i, j

    ! To each partial CF cell, assign the maximum value from a full edge-adjacent cell
    do j = 2, ny-1
       do i = 2, nx-1
          if (partial_cf_mask(i,j) == 1) then
             field1(i,j) = &
                  max(field1(i-1,j)*full_mask(i-1,j), field1(i+1,j)*full_mask(i+1,j), &
                      field1(i,j-1)*full_mask(i,j-1), field1(i,j+1)*full_mask(i,j+1))
             if (present(field2)) then
                field2(i,j) = &
                     max(field2(i-1,j)*full_mask(i-1,j), field2(i+1,j)*full_mask(i+1,j), &
                         field2(i,j-1)*full_mask(i,j-1), field2(i,j+1)*full_mask(i,j+1))
             endif
          endif
       enddo
    enddo

  end subroutine glissade_extrapolate_to_calving_front

!---------------------------------------------------------------------------

end module glissade_calving

!---------------------------------------------------------------------------
