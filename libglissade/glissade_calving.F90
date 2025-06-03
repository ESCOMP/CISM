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

  use glimmer_paramets, only: eps11, thk0
  use glimmer_physcon, only: rhoi, rhoo, grav, scyr
  use glide_diagnostics, only: point_diag

  implicit none

  private
  public :: glissade_calving_mask_init, glissade_calve_ice, &
            glissade_remove_icebergs, glissade_remove_isthmuses, glissade_limit_cliffs,  &
            glissade_stress_tensor_eigenvalues, glissade_strain_rate_tensor_eigenvalues
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

  subroutine glissade_calve_ice(nx,             ny,      &
                                which_calving,           &
                                calving_domain,          &
                                which_ho_calving_front,  &
                                which_ho_calvingmip_domain, &
                                parallel,                &
                                calving,                 &  ! calving derived type
                                itest,  jtest,  rtest,   &
                                dt,             time,    &  ! s
                                dx,             dy,      &  ! m
                                x0,             y0,      &  ! m
                                x1,             y1,      &  ! m
                                sigma,                   &
                                thklim,                  &  ! m
                                uvel_2d,        vvel_2d, &  ! m/s
                                thck_pre_transport,      &  ! m
                                thck,           relx,    &  ! m
                                topg,           eus)        ! m

    ! Calve ice according to one of several methods.
    ! Note: This subroutine uses SI units.

    !use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask
    use glissade_masks, only: glissade_get_masks, glissade_marine_cliff_mask, glissade_marine_grounded_mask, &
         glissade_grounded_calving_front_mask, glissade_secondary_floating_mask
    use glissade_utils, only: glissade_input_fluxes
    use glissade_grid_operators, only: glissade_unstagger

    implicit none

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

    integer, intent(in) :: nx, ny                  !> horizontal grid dimensions

    !TODO: Move these options to the calving derived type
    integer, intent(in) :: which_calving           !> option for calving law
    integer, intent(in) :: calving_domain          !> option for where calving can occur
                                                   !> = 0 if calving occurs at the ocean edge only
                                                   !> = 1 if calving occurs everywhere the calving criterion is met
                                                   !> = 2 if calving occurs where criterion is met and there is a connected path
                                                   !>     to the ocean through other cells where the criterion is met
    integer, intent(in) :: which_ho_calving_front  !> = 1 for subgrid calving-front scheme, else = 0
    integer, intent(in) :: which_ho_calvingmip_domain  !> = 1 for circular, 2 for Thule; otherwise = 0

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
!    real(dp), dimension(:,:), intent(inout)  :: effective_areafrac  !> effective fractional area, < 1 for partial CF cells
!    real(dp), dimension(:,:), intent(inout)  :: lateral_rate        !> lateral calving rate (m/s) at calving front
!    real(dp), dimension(:,:), intent(in)     :: tau_eigen1          !> first eigenvalue of 2D horizontal stress tensor (Pa)
!    real(dp), dimension(:,:), intent(in)     :: tau_eigen2          !> second eigenvalue of 2D horizontal stress tensor (Pa)
!    real(dp), dimension(:,:), intent(in)     :: eps_eigen1          !> first eigenvalue of 2D horizontal strain-rate tensor (1/s)
!    real(dp), dimension(:,:), intent(in)     :: eps_eigen2          !> second eigenvalue of 2D horizontal strain-rate tensor (1/s)
!    real(dp), dimension(:,:), intent(in)     :: thermal_forcing_applied,& ! thermal_forcing applied (degC)
!    real(dp), dimension(:,:), intent(in)     :: runoff_applied      ! runoff_applied (kg/s/m2) ISMIP6 standard, divide by 1000 to get m/s
!    real(dp), dimension(:,:,:), intent(inout):: damage              !> 3D scalar damage parameter
!    real(dp), intent(in)                     :: damage_threshold    !> threshold value where ice is sufficiently damaged to calve
!    real(dp), intent(in)                     :: damage_constant     !> rate of change of damage (1/s) per unit stress (Pa)
!    real(dp), intent(in)            :: cf_advance_retreat_amplitude !> amplitude (m/yr) of CF advance/retreat rate
!    real(dp), intent(in)            :: cf_advance_retreat_period    !> period (yr) of CF advance/retreat rate
!    integer,  dimension(:,:), intent(in)     :: protected_mask      !> integer mask: = 1 for cells that are able to fill,
!                                                                    !>  = 0 for cells that are not allowed to fill
!    integer,  dimension(:,:), intent(in)     :: damage_mask         !> integer mask: = 1 for damaged cells, else = 0
!    integer,  dimension(:,:), intent(in)     :: calving_mask        !> integer mask: calve ice where calving_mask = 1
!    real(dp), dimension(:,:), intent(out)    :: calving_thck        !> thickness lost due to calving in each grid cell (m)

    integer, intent(in) :: itest, jtest, rtest                     !> coordinates of diagnostic point
    real(dp), intent(in)                      :: dt                !> model timestep (s)
    real(dp), intent(in)                      :: time              !> model time (s)
    real(dp), intent(in)                      :: dx, dy            !> grid cell size in x and y directions (m)
    real(dp), dimension(nx-1), intent(in)     :: x0                !> x coordinates of NE cell corners (m)
    real(dp), dimension(ny-1), intent(in)     :: y0                !> y coordinates of NE cell corners (m)
    real(dp), dimension(nx), intent(in)       :: x1                !> x coordinates of cell centers (m)
    real(dp), dimension(ny), intent(in)       :: y1                !> y coordinates of cell centers (m)
    real(dp), dimension(:), intent(in)        :: sigma             !> vertical sigma coordinate
    real(dp), dimension(nx-1,ny-1), intent(in):: uvel_2d, vvel_2d  !> mean ice velocity components at vertices (m/s)
    real(dp), dimension(nx,ny), intent(in)    :: thck_pre_transport!> ice thickness (m) before doing transport, SMB, and BMB
    real(dp), dimension(nx,ny), intent(inout) :: thck              !> ice thickness (m)
    real(dp), dimension(nx,ny), intent(in)    :: relx              !> relaxed bedrock topography (m)
    real(dp), dimension(nx,ny), intent(in)    :: topg              !> present bedrock topography (m)
    real(dp), intent(in)                      :: thklim            !> minimum thickness for dynamically active grounded ice (m)
    real(dp), intent(in)                      :: eus               !> eustatic sea level (m)

    ! local variables

    integer :: nz          ! number of vertical levels
                           ! Note: number of ice layers = nz-1
    integer :: i, j, k, n, ig, jg
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
         marine_grounded_mask,   & ! = 1 where ice is grounded and marine-based 
         calving_front_mask,     & ! = 1 where ice is floating with at least one ocean edge neighbor, else = 0
         secondary_floating_mask   ! = 1 where ice is floating and not connected to marine grounded ice, else = 0

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
         calving_dthck,        & ! thickness increment (m) to be added to calving%thck
         cf_length               ! length of calving front within a cell

    integer, dimension(nx,ny) :: &
         partial_cf_mask,      & ! = 1 for partially filled CF cells (thck < thck_effective), else = 0
         full_mask               ! = 1 for ice-filled cells that are not partial_cf cells, else = 0

    real(dp), dimension(-1:1,-1:1,nx,ny) :: &
         flux_in                 ! ice volume fluxes (m^3/s) into cell from each neighbor cell

    real(dp), dimension(nx-1,ny-1) :: &
         velnorm_mean            !  mean ice speed at vertices (m/s)

    integer, dimension(nx-1,ny-1) :: &
         vmask                   ! = 1 for vertices of active cells

    real(dp), dimension(nx,ny) :: &
         speed                   ! 2D ice speed averaged to cell centers (m/s)

    real(dp), dimension(2,8) :: &
         cf_location             ! x and y components of calving front location
                                 ! first index is (x,y); second corresponds to 8 CalvingMIP axes

    ! some optional diagnostics
    real(dp) :: &
         total_ice_area,       & ! total effective ice area (with weighting by effective_areafrac)
         total_cf_length         ! total length of the calving front

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

    !WHL - Not sure if this update is needed
    call parallel_halo(thck, parallel)

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

    ! Combine CALVING_FLOAT_ZERO with other (subgrid CF) options, here CF_ADVANCE_RETREAT_RATE and similar
    !   The idea is to remove flaoting ice first and then apply calving processes to the marine grounded ice
    !   and, in the subgrid case, partial cells adjacent to those
    ! TODO: How does this work together with subgrid protected cells?
    !       Is it more agressive to float kill first, like done here?
    
    if (which_calving == CF_ARR_FLOAT_ZERO .or. which_calving == CF_MR_FLOAT_ZERO .or. which_calving == CF_SLATER_MR_FLOAT_ZERO) then

       call parallel_halo(thck, parallel)

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

       ! Only for diagnostic here
       call glissade_marine_cliff_mask(&
            nx,            ny,                &
            ice_mask,      floating_mask,     &
            land_mask,     ocean_mask,        &
            calving%marine_cliff_mask)
       
       call parallel_halo(calving%marine_cliff_mask, parallel)    
       
       ! set the calving-law mask
       ! Note: Cells that meet the calving-law criteria will be calved provided they also lie in the calving domain,
       !       as determined below.
   
       if (which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then
          
          ! In the subgrid case, remove secondary floating ice only when not protected
          ! Rely on protected cell mechanism       
          call glissade_marine_grounded_mask(&
               nx,            ny,             &
               ice_mask,      floating_mask,  &
               land_mask,                     &
               marine_grounded_mask)
   
          call parallel_halo(marine_grounded_mask, parallel)    
   
          call glissade_secondary_floating_mask(&
               nx,            ny,             &
               floating_mask,                 &
               marine_grounded_mask,          &
               secondary_floating_mask)
          
          call parallel_halo(secondary_floating_mask, parallel)    
   
          do j = 1, ny
             do i = 1, nx
                ! removing only floating ice that is not connected to a marine grounded cell,
                !  keeping cells that become subgrid calving front cells
                if (secondary_floating_mask(i,j) == 1) then
                   calving_law_mask(i,j) = .true.
                else
                   calving_law_mask(i,j) = .false.
                endif
             enddo
          enddo
   
       else ! other calving options (no subgrid calving front)

          do j = 1, ny
             do i = 1, nx
                ! Remove all floating ice
                if (floating_mask(i,j) == 1) then
                   calving_law_mask(i,j) = .true.
                else
                   calving_law_mask(i,j) = .false.
                endif
             enddo
          enddo
   
       endif ! which_ho_calving_front

       
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
!!                   write(6,*) 'Helo Calve ice: task, i, j, calving_thck =', this_rank, i, j, float_fraction_calve * thck(i,j)
                endif

                calving%calving_thck(i,j) = calving%calving_thck(i,j) + float_fraction_calve * thck(i,j)
                thck(i,j) = thck(i,j) - float_fraction_calve * thck(i,j)
            endif
          enddo
       enddo
       
    endif ! which_calving == CF_ARR_FLOAT_ZERO .or. which_calving == CF_MR_FLOAT_ZERO for combination
    
    ! On with physical calving options
    
    ! Calving schemes with a subgrid calving front:
    
    if (which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then

       !TODO - Rename 'protected', since unprotected cells are in fact protected now.
       !
       ! Use one of the subgrid calving front schemes:
       ! * prescribed advance/retreat rate
       ! * thickness-based calving
       ! * eigencalving
       ! * damage-based calving
       ! Each of these follows a similar pattern:
       ! (1) Where ice has been transported downstream from a partial CF cell
       !     to previously ice-free cells, move it back upstream.
       ! (2) Compute some masks related to calving.
       ! (3) Depending on the calving law, compute the lateral calving rate
       !     and convert to a thinning rate.
       ! (4) Apply the calving-derived thinning. If a full column is removed at the CF,
       !     do additional thinning upstream.
       ! (5) Where H > H_effective in CF cells, set H = H_effective and
       !     move the extra ice downstream, advancing the CF.
       ! Only step (3) depends on the specific calving law.
       
       ! Compute the ice speed at cell centers, averaged from neighboring vertices.
       ! Include in the average only vertices with nonzero speeds (i.e., ice present)
       ! This speed is used to compute the calving rate for thickness-based calving,
       !  eigencalving, and damage-based calving.

       velnorm_mean = sqrt(uvel_2d**2 + vvel_2d**2)

       where (velnorm_mean > 0.0d0)
          vmask = 1
       elsewhere
          vmask = 0
       endwhere

       ! Interpolate the speed from cell vertices to centers.
       ! 'stagger_margin_in = 1' means that masked-out values are not part of the average.

       call glissade_unstagger(&
            nx,            ny,      &
            velnorm_mean,  speed,   &
            vmask,         stagger_margin_in = 1)

       call parallel_halo(speed, parallel)

       ! Compute the ice flux into each cell from each neighbor cell.
       ! This is an upwind estimate based on cell-center thickness.
       ! It is not equivalent to computing the incremental remapping flux,
       !  but near the ice edge (where reconstructed thicknesses near cell edges
       !  are close to cell-center values) it is a good approximation.

       call glissade_input_fluxes(&
            nx,      ny,                      &
            dx,      dy,                      &
            itest,   jtest,  rtest,           &
            thck_pre_transport,               & ! m
            uvel_2d, vvel_2d,                 & ! m/s
            flux_in,                          & ! m^3/s
            parallel)

       if (verbose_calving) then
          call point_diag(thck, 'Before redistribution, thck (m)', itest, jtest, rtest, 7, 7)
       endif

       ! Gather ice that has flowed to unprotected cells and move it back upstream

       call redistribute_unprotected_ice(&
            nx,                ny,            &
            itest,  jtest,  rtest,            &
            parallel,                         &
            calving%protected_mask,           &
            flux_in,                          & ! m^3/s
            thck,                             & ! m
            calving%calving_thck)               ! m

       if (verbose_calving) then
          call point_diag(thck, 'After redistribution, thck (m)', itest, jtest, rtest, 7, 7)
       endif

       ! Compute masks for calving.

       call glissade_get_masks(&
            nx,            ny,             &
            parallel,                      &
            thck,          topg,           &
            eus,           thklim,         &
            ice_mask,                      &
            floating_mask = floating_mask, &
            ocean_mask = ocean_mask,       &
            land_mask = land_mask)

! commented to use the grounded calving front instead       
!       call glissade_calving_front_mask(&
!            nx,            ny,             &
!            which_ho_calving_front,        &
!            parallel,                      &
!            thck,          topg,           &
!            eus,                           &
!            ice_mask,      floating_mask,  &
!            ocean_mask,    land_mask,      &
!            calving_front_mask,            &
!            calving%dthck_dx_cf,           &
!            dx,            dy,             &
!            calving%thck_effective,        &
!            calving%thck_effective_min,    &
!            partial_cf_mask,               &
!            full_mask,                     &
!            calving%effective_areafrac)

       call glissade_grounded_calving_front_mask(&
            nx,                     ny,    &
            which_ho_calving_front,        &
            parallel,                      &
            thck,                   topg,  &
            eus,                           &
            ice_mask,      floating_mask,  &
            ocean_mask,        land_mask,  &
            calving_front_mask,            &
            calving%dthck_dx_cf,           &
            dx,            dy,             &
            calving%thck_effective,        &
            calving%thck_effective_min,    &
            partial_cf_mask,               &
            full_mask,                     &
            calving%effective_areafrac)

       if (verbose_calving) then
          call point_diag(calving_front_mask, 'calving_front_mask', itest, jtest, rtest, 7, 7)
          call point_diag(partial_cf_mask, 'partial_cf_mask', itest, jtest, rtest, 7, 7)
          call point_diag(full_mask, 'full_mask', itest, jtest, rtest, 7, 7)
          call point_diag(calving%thck_effective, 'thck_effective (m)', itest, jtest, rtest, 7, 7)
          call point_diag(calving%effective_areafrac, 'effective_areafrac', itest, jtest, rtest, 7, 7)
       endif

       ! Compute the effective length of the calving front in each grid cell

       if (which_calving == CF_ADVANCE_RETREAT_RATE) then

          ! compute the CF length as a function of a cell's location on the unit circle
          ! surrounding the origin (assuming a radially symmetric calving rate).

          call compute_calving_front_length_radial(&
               nx,           ny,             &
               dx,           dy,             &
               x1,           y1,             &
               itest, jtest, rtest,          &
               calving_front_mask,           &
               ocean_mask,                   &
               cf_length)

       else

          ! compute the CF length for each cell based on its number of ocean neighbors
          ! (i.e., enhanced calving for cells with 2 or 3 ocean neighbors).

          call compute_calving_front_length(&
               nx,           ny,             &
               dx,           dy,             &
               itest, jtest, rtest,          &
               calving_front_mask,           &
               ocean_mask,                   &
               cf_length)

       endif   ! which_calving

       call parallel_halo(cf_length, parallel)

       if (verbose_calving) then
          call point_diag(cf_length, 'cf_length (m)', itest, jtest, rtest, 7, 7)
          ! Diagnose the total CF length
          total_cf_length = 0.0d0
          do j = nhalo+1, ny-nhalo
             do i = nhalo+1, nx-nhalo
                if (calving_front_mask(i,j) == 1) then
                   total_cf_length = total_cf_length + cf_length(i,j)
                endif
             enddo
          enddo
          total_cf_length =  parallel_reduce_sum(total_cf_length)
          if (this_rank == rtest) then
             print*, 'Total CF length (km)', total_cf_length/1000.d0
          endif
       endif

       ! Depending on the calving method, compute the calving rate for each grid cell
       ! and convert to an equivalent thinning rate.

!HG<       if (which_calving == CF_ADVANCE_RETREAT_RATE) then
       if (which_calving == CF_ADVANCE_RETREAT_RATE .or. which_calving == CF_ARR_FLOAT_ZERO) then

          call calving_front_advance_retreat(&
               nx,                 ny,                    &
               dx,                 dy,                    &
               dt,                 time,                  &  ! s
               itest,   jtest,     rtest,                 &
               calving_front_mask,                        &
               thck_pre_transport,                        &  ! m
               thck,                                      &  ! m
               cf_length,                                 &  ! m
               calving%thck_effective,                    &  ! m
               calving%cf_advance_retreat_amplitude/scyr, &  ! m/s
               calving%cf_advance_retreat_period*scyr,    &  ! s
               calving_dthck)                                ! m

       elseif (which_calving == CF_MELTRATE .or. which_calving == CF_MR_FLOAT_ZERO) then

          call calving_front_meltrate(&
               nx,                 ny,                    &
               dx,                 dy,                    &
               dt,                 time,                  &  ! s
               itest,   jtest,     rtest,                 &
               calving_front_mask,                        &
               thck_pre_transport,                        &  ! m
               thck,                                      &  ! m
               cf_length,                                 &  ! m
               calving%thck_effective,                    &  ! m
               calving%cf_advance_retreat_amplitude/scyr, &  ! m/s
               calving%cf_advance_retreat_period*scyr,    &  ! s
               calving_dthck)                                ! m

       elseif (which_calving == CF_SLATER) then
          
          call slater_calving_meltrate(&
               nx,                 ny,                    &
               dx,                 dy,                    &
               dt,                 time,                  &  ! s
               itest,   jtest,     rtest,                 &
               calving_front_mask,                        &
               thck_pre_transport,                        &  ! m
               calving%thermal_forcing_applied,           &  ! degC
               calving%runoff_applied/1000.,              &  ! m/s; (read in kg/m2/s) 
               thck,                                      &  ! m
               cf_length,                                 &  ! m
               calving%thck_effective,                    &  ! m
               calving%cf_advance_retreat_amplitude/scyr, &  ! m/s
               calving%cf_advance_retreat_period*scyr,    &  ! s
               calving_dthck)                                ! m

       elseif (which_calving == CF_SLATER_MR .or. which_calving == CF_SLATER_MR_FLOAT_ZERO) then

          call slater_meltrate(&
               nx,                 ny,                    &
               dx,                 dy,                    &
               dt,                 time,                  &  ! s
               itest,   jtest,     rtest,                 &
               calving_front_mask,                        &
               thck_pre_transport,                        &  ! m
               calving%thermal_forcing_applied,           &  ! degC
               calving%runoff_applied/1000.,              &  ! m/s; (read in kg/m2/s) 
               thck,                                      &  ! m
               cf_length,                                 &  ! m
               calving%thck_effective,                    &  ! m
               calving%cf_advance_retreat_amplitude/scyr, &  ! m/s
               calving%cf_advance_retreat_period*scyr,    &  ! s
               calving_dthck)                                ! m

       elseif (which_calving == CALVING_THCK_THRESHOLD) then

          call extrapolate_to_calving_front(&
               nx,                 ny,     &
               partial_cf_mask,            &
               full_mask,                  &
               calving%effective_areafrac, &
               speed)

          call thickness_based_calving(&
               nx,                 ny,                    &
               dx,                 dy,                    &  ! m
               dt,                                        &  ! s
               itest,   jtest,     rtest,                 &
               calving_front_mask,                        &
               speed,                                     &  ! m/s
               cf_length,                                 &  ! m
               calving%thck_effective,                    &  ! m
               calving%minthck,                           &  ! m
               calving_dthck)                                ! m

       elseif (which_calving == CALVING_STRESS) then

          call extrapolate_to_calving_front(&
               nx,                 ny,     &
               partial_cf_mask,            &
               full_mask,                  &
               calving%effective_areafrac, &
               speed)

          call extrapolate_to_calving_front(&
               nx,                 ny,     &
               partial_cf_mask,            &
               full_mask,                  &
               calving%effective_areafrac, &
               calving%tau_eigen1,         &
               calving%tau_eigen2)

          call stress_based_calving(&
               nx,                 ny,            &
               dx,                 dy,            &  ! m
               dt,                                &  ! s
               itest,   jtest,     rtest,         &
               calving_front_mask,                &
               speed,                             &  ! m/s
               cf_length,                         &  ! m
               calving%thck_effective,            &  ! m
               calving%tau_eigen1,                &  ! Pa
               calving%tau_eigen2,                &  ! Pa
               calving%tau_eigenconstant1,        &
               calving%tau_eigenconstant2,        &
               calving%stress_threshold,          &  ! Pa
               calving%lateral_rate_min/scyr,     &  ! m/s
               calving_dthck)                        ! m

       elseif (which_calving == CALVING_STRESS_STOCHASTIC) then

          call extrapolate_to_calving_front(&
               nx,                 ny,     &
               partial_cf_mask,            &
               full_mask,                  &
               calving%effective_areafrac, &
               speed)

          call extrapolate_to_calving_front(&
               nx,                 ny,     &
               partial_cf_mask,            &
               full_mask,                  &
               calving%effective_areafrac, &
               calving%tau_eigen1,         &
               calving%tau_eigen2)

          call stochastic_stress_based_calving(&
               nx,                 ny,            &
               dx,                 dy,            &  ! m
               dt,                                &  ! s
               itest,   jtest,     rtest,         &
               parallel,                          &
               calving_front_mask,                &
               thck,                              &  ! m
               calving%thck_effective,            &  ! m
               calving%effective_areafrac,        &
               speed,                             &  ! m/s
               cf_length,                         &  ! m
               calving%tau_eigen1,                &  ! Pa
               calving%tau_eigen2,                &  ! Pa
               calving%tau_eigenconstant1,        &
               calving%tau_eigenconstant2,        &
               calving%stress_threshold,          &  ! Pa
               calving%effec_stress_min,          &  ! Pa
               calving%length_scale,              &  ! m
               calving_dthck)                        ! m

       elseif (which_calving == EIGEN_CALVING) then

          call extrapolate_to_calving_front(&
               nx,                 ny,     &
               partial_cf_mask,            &
               full_mask,                  &
               calving%effective_areafrac, &
               calving%eps_eigen1,         &
               calving%eps_eigen2)

          call eigencalving(&
               nx,                 ny,            &
               dx,                 dy,            &  ! m
               dt,                                &  ! s
               itest,   jtest,     rtest,         &
               calving_front_mask,                &
               cf_length,                         &  ! m
               calving%thck_effective,            &  ! m
               calving%eps_eigen1,                &  ! 1/s
               calving%eps_eigen2,                &  ! 1/s
               calving%eigenconstant,             &  ! m
               calving_dthck)                        ! m

       elseif (which_calving == CALVING_DAMAGE) then

          call extrapolate_to_calving_front(&
               nx,                 ny,     &
               partial_cf_mask,            &
               full_mask,                  &
               calving%effective_areafrac, &
               speed)

          !Note - Optionally, instead of weighting by areafrac, we could assign zero weight
          !        to the CF cell and use the upstream value. This would increase damage at the CF.

          call extrapolate_to_calving_front(&
               nx,                 ny,     &
               partial_cf_mask,            &
               full_mask,                  &
               calving%effective_areafrac, &
               calving%tau_eigen1,         &
               calving%tau_eigen2)

          call stochastic_damage_based_calving(&
               nx,       ny,       nz,                 &
               dx,                 dy,                 &  ! m
               sigma,              dt,                 &
               itest,   jtest,     rtest,              &
               parallel,                               &
               floating_mask,                          &
               calving_front_mask,                     &
               thck,                                   &  ! m
               topg,                                   &  ! m
               calving%tau_eigen1, calving%tau_eigen2, &  ! Pa
               calving%tau_eigenconstant1,             &
               calving%tau_eigenconstant2,             &
               calving%stress_threshold,               &
               calving%effec_stress_min,               &  ! Pa
               calving%damage_constant*scyr,           &  ! Pa s
               calving%damage,                         &
               calving_dthck,                          &  ! m
               calving%eps_eigen1, calving%eps_eigen2)    ! 1/s

       endif   ! which_calving

       call parallel_halo(calving_dthck, parallel)

       ! Apply calving_dthck as computed above.

       if (which_calving == CALVING_DAMAGE) then

          ! different treatment because we can calve cells not on the CF
          where (calving_dthck == thck)
             calving%calving_thck = calving_dthck
             thck = 0.0d0
          endwhere

       else

!          call apply_calving_dthck(&
!               nx,           ny,        &
!               itest, jtest, rtest,     &
!               parallel,                &
!               calving_front_mask,      &
!               floating_mask,           &
!               full_mask,               &
!               flux_in,                 &
!               calving_dthck,           &
!               thck,                    &
!               calving%calving_thck)

          call apply_grounded_calving_dthck(&
               nx,           ny,        &
               itest, jtest, rtest,     &
               parallel,                &
               calving_front_mask, &
               ice_mask,                &
               floating_mask,           &
               full_mask,               &
               flux_in,                 &
               calving_dthck,           &
               thck,                    &
               calving%calving_thck)

       endif

       !TODO - Add a bug check for negative thicknesses?
       thck = max(thck, 0.0d0)

       call parallel_halo(thck, parallel)
       call parallel_halo(calving%calving_thck, parallel)

       ! Recompute the calving masks

       call glissade_get_masks(&
            nx,            ny,             &
            parallel,                      &
            thck,          topg,           &
            eus,           thklim,         &
            ice_mask,                      &
            floating_mask = floating_mask, &
            ocean_mask = ocean_mask,       &
            land_mask = land_mask)

!       call glissade_calving_front_mask(&
!            nx,            ny,             &
!            which_ho_calving_front,        &
!            parallel,                      &
!            thck,          topg,           &
!            eus,                           &
!            ice_mask,      floating_mask,  &
!            ocean_mask,    land_mask,      &
!            calving_front_mask,            &
!            calving%dthck_dx_cf,           &
!            dx,            dy,             &
!            calving%thck_effective,        &
!            calving%thck_effective_min,    &
!            partial_cf_mask,               &
!            full_mask,                     &
!            calving%effective_areafrac)

       call glissade_grounded_calving_front_mask(&
            nx,                     ny,    &
            which_ho_calving_front,        &
            parallel,                      &
            thck,                   topg,  &
            eus,                           &
            ice_mask,      floating_mask,  &
            ocean_mask,        land_mask,  &
            calving_front_mask,            &
            calving%dthck_dx_cf,           &
            dx,            dy,             &
            calving%thck_effective,        &
            calving%thck_effective_min,    &
            partial_cf_mask,               &
            full_mask,                     &
            calving%effective_areafrac)

       if (which_calving == CF_ADVANCE_RETREAT_RATE) then

          ! Compute some CalvingMIP diagnostics.
          ! Note: If running with a prescribed advance/retreat rate on a grid other than
          !       the CalvingMIP circular and Thule domains, we would need some additional
          !       logic to identify the domain.
          ! Note: With a prescribed advance/retreat rate, it isn't necessary to call advance_calving_front.

          if (which_ho_calvingmip_domain == HO_CALVINGMIP_DOMAIN_CIRCULAR) then

             call locate_calving_front_circular(&
                  nx,             ny,           &
                  dx,             dy,           &  ! m
                  x0,             y0,           &  ! m
                  x1,             y1,           &  ! m
                  parallel,                     &
                  itest, jtest, rtest,          &
                  calving%effective_areafrac,   &
                  cf_location)                     ! m

          elseif (which_ho_calvingmip_domain == HO_CALVINGMIP_DOMAIN_THULE) then

                call locate_calving_front_thule(&
                  nx,             ny,           &
                  dx,             dy,           &  ! m
                  x0,             y0,           &  ! m
                  x1,             y1,           &  ! m
                  parallel,                     &
                  itest, jtest, rtest,          &
                  calving%effective_areafrac,   &
                  cf_location)                     ! m

          endif

          ! Compute the total ice area and the area of each quadrant
          total_ice_area = 0.0d0
          do j = nhalo+1, ny-nhalo
             do i = nhalo+1, nx-nhalo
                total_ice_area = total_ice_area + dx*dy*calving%effective_areafrac(i,j)
             enddo
          enddo
          total_ice_area = parallel_reduce_sum(total_ice_area)
          if (this_rank == rtest) then
             print*, 'Total ice area (km^2)=', total_ice_area/1.0d6
             print*, 'Quadrant area (km^2)=', total_ice_area/4.0d6
          endif

          if (verbose_calving) then
             call point_diag(calving%thck_effective, 'New thck_effective (m)', itest, jtest, rtest, 7, 7)
             call point_diag(calving%effective_areafrac, 'New effective_areafrac', itest, jtest, rtest, 7, 7)
!!             call point_diag(thck/calving%effective_areafrac, '   CF thickness', itest, jtest, rtest, 7, 7)
          endif

       else  ! other calving schemes

          ! Where thck > thck_effective, allow the CF to advance by distributing ice downstream.
          ! Note: This is not necessary when running with a prescribed advance/retreat rate.

          call advance_calving_front(&
               nx,           ny,        &
               itest, jtest, rtest,     &
               ocean_mask,              &
               calving_front_mask,      &
               flux_in,                 &
               calving%thck_effective,  &
               thck)

          ! Recompute the calving masks
          !TODO - Are these calls needed? These subroutines are called again before the velocity solver.

          call glissade_get_masks(&
               nx,            ny,             &
               parallel,                      &
               thck,          topg,           &
               eus,           thklim,         &
               ice_mask,                      &
               floating_mask = floating_mask, &
               ocean_mask = ocean_mask,       &
               land_mask = land_mask)

!          call glissade_calving_front_mask(&
!               nx,            ny,             &
!               which_ho_calving_front,        &
!               parallel,                      &
!               thck,          topg,           &
!               eus,                           &
!               ice_mask,      floating_mask,  &
!               ocean_mask,    land_mask,      &
!               calving_front_mask,            &
!               calving%dthck_dx_cf,           &
!               dx,            dy,             &
!               calving%thck_effective,        &
!               calving%thck_effective_min,    &
!               partial_cf_mask,               &
!               full_mask,                     &
!               calving%effective_areafrac)

          call glissade_grounded_calving_front_mask(&
               nx,                     ny,    &
               which_ho_calving_front,        &
               parallel,                      &
               thck,                   topg,  &
               eus,                           &
               ice_mask,      floating_mask,  &
               ocean_mask,        land_mask,  &
               calving_front_mask,            &
               calving%dthck_dx_cf,           &
               dx,            dy,             &
               calving%thck_effective,        &
               calving%thck_effective_min,    &
               partial_cf_mask,               &
               full_mask,                     &
               calving%effective_areafrac)
       
       endif

    else   ! other calving options (no subgrid calving front)
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

    endif   ! which_ho_calving_front


    if (verbose_calving) then
       call point_diag(thck, 'After calving, new thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(calving%calving_thck, 'calving_thck (m)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine glissade_calve_ice

!---------------------------------------------------------------------------

  subroutine redistribute_unprotected_ice(&
       nx,              ny,       &
       itest,  jtest,   rtest,    &
       parallel,                  &
       protected_mask,            &
       flux_in,                   &
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

    real(dp), dimension(-1:1,-1:1,nx,ny), intent(in) :: &
         flux_in                  ! ice volume fluxes (m^3/s) into cell from each neighbor cell

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m) before and after redistribution
         calving_thck              ! thickness (m) calved from each cell

    ! local variables

    integer :: i, j, ii, jj, count
    integer :: iup, jup, idn, jdn
    real(dp) :: thck_max           ! max thickness (m) in protected upstream cell
    real(dp) :: dthck              ! ice thickness (m) to be redistributed
    real(dp) :: total_flux         ! total flux (m^3/s) entering a cell from neighbor cells
    real(dp) :: total_dthck        ! total thickness (m) to be redistributed

    ! Write the initial thicknesses
    ! I think the halo update is needed only to get the right halo values for diagnostics
    call parallel_halo(thck, parallel)

    ! Identify unprotected ice with nonzero thickness.
    ! Instead of calving this ice, move it to one or more protected upstream CF cell
    ! (from which most or all of the ice likely arrived during transport).

    do j = 2, ny-1
       do i = 2, nx-1
          if (thck(i,j) > 0.0d0 .and. protected_mask(i,j) == 0) then

             ! Given flux_in (ice flux in m^3/s entering the cell from each upstream CF neighbor),
             ! compute the fraction of the flux to give back to each upstream neighbor.
             count = 0
             total_flux = 0.0d0
             do jj = -1, 1
                do ii = -1, 1
                   iup = i + ii; jup = j + jj
                   if (flux_in(ii,jj,i,j) > 0.0d0 .and. protected_mask(iup,jup) == 1) then
                      count = count + 1
                      total_flux = total_flux + flux_in(ii,jj,i,j)
                   endif
                enddo
             enddo

             ! Move ice from the unprotected cell to its upstream neighbors.
             ! This can result in H > H_eff in upstream cells, but the excess ice will be removed
             !  later by calving or downstream redistribution.
             total_dthck = thck(i,j)
             if (total_flux > 0.0d0) then
                do jj = -1, 1
                   do ii = -1, 1
                      iup = i + ii; jup = j + jj
                      if (flux_in(ii,jj,i,j) > 0.0d0 .and. protected_mask(iup,jup) == 1) then
                         dthck = total_dthck * flux_in(ii,jj,i,j)/total_flux
                         thck(iup,jup) = thck(iup,jup) + dthck
                         thck(i,j) = thck(i,j) - dthck
!                         if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
!                            print*, '   Upstream ii, jj, frac:', ii, jj, flux_in(ii,jj,i,j)/total_flux
!                         endif
                      endif
                   enddo
                enddo
             endif

          endif   ! thck > 0 and unprotected
       enddo   ! i
    enddo   ! j

    call parallel_halo(thck, parallel)

  end subroutine redistribute_unprotected_ice

!---------------------------------------------------------------------------

  subroutine compute_calving_front_length(&
       nx,           ny,             &
       dx,           dy,             &
       itest, jtest, rtest,          &
       calving_front_mask,           &
       ocean_mask,                   &
       cf_length)

    ! Compute the effective length of the calving front in each grid cell,
    !  based on the number of ocean neighbors.
    ! Cells with a single ocean edge receive a length of dx or dy.
    ! Cells with two or three adjacent edges receive a length of sqrt(dx^2 + dy^2).
    !
    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy                    ! grid cell size (m)

    integer, dimension(nx,ny), intent(in) :: &
         calving_front_mask,     & ! = 1 for floating cells with an ice-free ocean neighber
         ocean_mask                ! = 1 for ice-free ocean cells

    real(dp), dimension(nx,ny), intent(out) :: &
         cf_length                 ! calving front length (m) through a given grid cell

    ! local variables

    integer :: i, j, ip, jp, ii, jj
    integer :: count, count_ew, count_ns    ! number of neighbors of a CF cell

    ! Compute the CF length based on the number of ocean neighbors;
    ! the CF is longer for cells with two ocean neighbors.

    cf_length = 0.0d0

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (calving_front_mask(i,j) == 1) then

             ! Count the number of ocean edges in each direction
             count_ns = ocean_mask(i,j-1) + ocean_mask(i,j+1)   ! N/S neighbors; edge of length dx
             count_ew = ocean_mask(i-1,j) + ocean_mask(i+1,j)   ! E/W neighbors; edge of length dy
             count = count_ns + count_ew
             if (count == 1) then   ! one ocean neighbor
                cf_length(i,j) = count_ns*dx + count_ew*dy
             elseif (count == 2) then   ! two ocean neighbors; usually a corner in the CF
                if (count_ew == 1 .and. count_ns == 1) then
                   cf_length(i,j) = sqrt(dx*dx + dy*dy)
                else  ! rare case of an isthmus with ocean on opposite edges
                   cf_length(i,j) = count_ns*dx + count_ew*dy
                endif
             elseif (count == 3) then
                ! This is a bit arbitrary; assume the same length as a CF bordering two edges.
                ! The goal is to speed up calving for 'teeth' but not remove them instantly.
                cf_length(i,j) = dx + dy
             endif

          endif   ! calving_mask
       enddo   ! i
    enddo   ! j

  end subroutine compute_calving_front_length

!---------------------------------------------------------------------------

  subroutine compute_calving_front_length_radial(&
       nx,           ny,             &
       dx,           dy,             &
       x1,           y1,             &
       itest, jtest, rtest,          &
       calving_front_mask,           &
       ocean_mask,                   &
       cf_length)

    ! This is an alternate method of computing the calving front length in idealized problems
    ! (e.g., CalvingMIP) in which the calving rate has radial symmetry.
    ! For this method, the CF length for a given gridcell has a minimum value = dx) when a line drawn
    !  from the origin to the cell center is parallel to the x or y axis.
    ! The CF length has a maximum value when a line drawn from the origin to the cell center
    !  makes an angle of pi/4 = 45 degrees with the x and y axes.
    ! At intermediate angles, the CF length has an intermediate value.

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy                    ! grid cell size (m)

    real(dp), dimension(nx), intent(in) :: x1  ! x coordinate of cell centers
    real(dp), dimension(ny), intent(in) :: y1  ! y coordinate of cell centers

    integer, dimension(nx,ny), intent(in) :: &
         calving_front_mask,     & ! = 1 for floating cells with an ice-free ocean neighber
         ocean_mask                ! = 1 for ice-free ocean cells

    real(dp), dimension(nx,ny), intent(out) :: &
         cf_length                 ! calving front length (m) through a given grid cell

    ! local variables

    integer :: i, j

    real(dp) :: absx, absy          ! absolute value of (x1,y1)
    real(dp) :: theta               ! angle relative to x or y axis

    if (abs(dx - dy) > eps11) then
       write(6,*) 'Error, compute_calving_front_length_radial, this_rank, i, j,', this_rank, i, j
       call write_log('Must have dx = dy for the cf_length method', GM_FATAL)
    endif

    ! Compute the CF length for each grid cell:
    ! *  equal to dx when the calving direction is parallel to the x or y axis
    ! *  equal to sqrt(2)*dx when the calving direction is along a diagonal
    ! *  of intermediate length for intermediate angles
    ! Note: The method assumes dx = dy.
    !       It fails if the center of a CF cell lies at the origin (0,0).

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (calving_front_mask(i,j) == 1) then
             absx = abs(x1(i))
             absy = abs(y1(j))
             if (absx == 0.0d0 .and. absy == 0.0d0) then
                write(6,*) 'Error, compute_calving_front_length_radial, this_rank, i, j,', this_rank, i, j
                call write_log('Cannot compute angle for CF cell at the origin', GM_FATAL)
             else
                if (absx >= absy) then
                   theta = atan(absy/absx)
                else   ! absx < absy
                   theta = atan(absx/absy)
                endif
             endif
             ! Note: theta lies in the range [0, pi/4], so cos(theta) is in the range [sqrt(2)/2, 1]
             cf_length(i,j) = dx / cos(theta)
          endif   ! calving front cell
       enddo   ! i
    enddo   ! j

  end subroutine compute_calving_front_length_radial

!---------------------------------------------------------------------------

  subroutine compute_calving_front_length_advance(&
       nx,           ny,             &
       dx,           dy,             &
       itest, jtest, rtest,          &
       calving_front_mask,           &
       cf_length)

    ! Compute the calving front length by a method slightly different from the method
    !  in subroutine compute_calving_front_length.
    ! This method is based on connecting all the adjacent cell centers on the calving front,
    !  such that cells with diagonal CF neighbors thin more than cells with edge neighbors.
    !
    ! This was originally the method used for CalvingMIP experiments,
    !  before switching to the more accurate radial method.
    ! When applied to CalvingMIP problems, it is more accurate than the method
    !  in subroutine compute_calving_front_length during the advance phases,
    !  but less accurate during the retreat phase.
    !
    ! A picture can be helpful. Let 'x' denote ice-filled cells and let 'o' denote ocean cells.
    !
    !  |    o   |    o   |    o   |    o   |
    !  |--------|--------|--------|--------|
    !  |        |        |        |        |
    !  |    x   |    x   |    o   |    o   |
    !  |        | (i,j+1)|        |        |
    !  |--------|--------|--------|--------|
    !  |        |        |        |        |
    !  |    x   |    x   |    x   |    x   |
    !  |        |  (i,j) | (i+1,j)|        |
    !  |------- |--------|--------|--------|
    !
    ! Note the corner in the CF at the NE vertex of cell (i,j).
    ! In general, we want calving to smooth out the CF rather than created teeth and indentations.
    ! During retreat, we can smooth the corner by enhancing the calving in cell (i,j+1),
    !  which borders two ocean cells. We do this by giving this cell an increased CF length.
    ! But when the CF is prescribed to advance, it's better to smooth the corner by reducing the
    !  amount of ice entering ocean cell (i+1,j+1). We do this by increasing the CF length in both
    !  (i,j+1) and (i+1,j), to reflect the diagonal connection between the two cells.
    ! This subroutine is not currently called but is here for now in case it's useful later.
    !---------------------------------------------------------------------------

    ! input/output variables

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy                    ! grid cell size (m)

    integer, dimension(nx,ny), intent(in) :: &
         calving_front_mask        ! = 1 for floating cells with an ice-free ocean neighber

    real(dp), dimension(nx,ny), intent(out) :: &
         cf_length                 ! calving front length (m) through a given grid cell

    ! local variables

    integer :: i, j, ip, jp, ii, jj
    integer :: count

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (calving_front_mask(i,j) == 1) then

             count = 0
             ! start with edge CF neighbors
             do jp = -1,1
                do ip = -1,1
                   ii = i + ip; jj = j + jp
                   if (calving_front_mask(ii,jj) == 1) then
                      if (abs(ip) == 1 .and. abs(jp) == 0) then ! E or W neighbor
                         count = count + 1
                         if (count <= 2) cf_length(i,j) = cf_length(i,j) + 0.5d0*dx
                      elseif (abs(jp) == 1 .and. abs(ip) == 0) then ! N or S neighbor
                         count = count + 1
                         if (count <= 2) cf_length(i,j) = cf_length(i,j) + 0.5d0*dy
                      endif
                   endif
                enddo   ! ip
             enddo   ! jp

             ! check for corner CF neighbors
             ! Note: Occasionally, a cell can have 3 CF neighbors.
             !      If so, then count corner neighbors only if there are not already 2 edge neighbors
             if (count < 2) then
                do jp = -1,1,2
                   do ip = -1,1,2
                      ii = i + ip; jj = j + jp
                      if (calving_front_mask(ii,jj) == 1) then
                         count = count + 1
                         if (count <= 2) then
                            cf_length(i,j) = cf_length(i,j) + 0.5d0*sqrt(dx*dx + dy*dy)
                         endif
                      endif
                   enddo
                enddo
             endif

             ! if we have not found 2 neighbors, then add a CF length of sqrt(dx*dy)/2 for each missing neighbor
             if (count == 0) then
                cf_length(i,j) = sqrt(dx*dy)
             elseif (count == 1) then
                cf_length(i,j) = cf_length(i,j) + 0.5d0*sqrt(dx*dy)
             endif

          endif   ! calving_mask
       enddo   ! i
    enddo   ! j

  end subroutine compute_calving_front_length_advance

!---------------------------------------------------------------------------

  subroutine thickness_based_calving(&
       nx,                 ny,                    &
       dx,                 dy,                    &  ! m
       dt,                                        &  ! s
       itest,   jtest,     rtest,                 &
       calving_front_mask,                        &
       speed,                                     &  ! m/s
       cf_length,                                 &  ! m
       thck_effective,                            &  ! m
       calving_minthck,                           &  ! m
       calving_dthck)                                ! m

    ! Calve ice that is thinner than a prescribed threshold.
    ! This option requires a subgrid calving front scheme.
    ! The CF thinning rate is based not on the nominal ice thickness (H = thck),
    !  but on the effective calving thickness (thck_effective = H_eff),
    !  which depends on the max thickness of the edge-adjacent floating cells.

    use glide_diagnostics, only: point_diag

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt                        ! time step (s)

    integer, dimension(nx,ny), intent(in)  ::  &
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         speed,                  & ! mean ice speed averaged to cell centers (m/s)
         cf_length,              & ! length of calving front in each grid cell (m)
         thck_effective            ! effective thickness for calving (m)

    real(dp), intent(in) :: &
         calving_minthck           ! target effective thickness (m) for the CF

    real(dp), dimension(nx,ny), intent(out) :: &
         calving_dthck             ! thickness lost due to calving (m)

    ! local variables

    integer :: i, j

    real(dp), dimension(nx,ny) :: &
         lateral_rate               ! lateral calving rate (m/s)

    ! Compute thinning in calving-front cells whose effective thickness (H_eff = thck_effective)
    !  is less than a prescribed minimum value (H_c = calving_minthck).
    !
    ! The lateral calving rate Cr is given by
    !
    !   Cr = max(0, 1 + (H_c - H_eff)/H_c) × v
    !
    ! where v is the ice speed (>=0) at the calving front.
    ! 
    ! In other words:
    ! * Cr = 0 where H_eff > 2*H_c
    ! * For H_eff < H < 2*H_c, there is some calving, but not enough to halt ice advance.
    ! * For H_eff = H_c, the calving rate is in balance with the ice speed, so the CF is stable.
    ! * For H_eff < H_c, the CF retreats. As H_eff -> 0, we have Cr -> 2v.
    !
    ! This calving law is based on Experiment 5 of CalvingMIP:
    !   https://github.com/JRowanJordan/CalvingMIP/wiki/Experiment-5
    !
    ! The lateral calving rate is converted to a thinning rate dH using
    ! 
    !   dH = Cr * dt * H_eff * cf_length / (dx*dy),
    !
    ! The RHS is equal to the ice volume removed per unit grid cell area during one time step

    lateral_rate = 0.0d0
    calving_dthck = 0.0d0

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (calving_front_mask(i,j) == 1) then

             lateral_rate(i,j) = &
                  max(0.0d0, 1.0d0 + (calving_minthck - thck_effective(i,j))/calving_minthck) * speed(i,j)
             calving_dthck(i,j) = lateral_rate(i,j) * dt * thck_effective(i,j) * cf_length(i,j) / (dx*dy)

          endif   ! CF mask
       enddo   ! i
    enddo   ! j

    if (verbose_calving) then
       call point_diag(speed*scyr, 'Thickness-based calving, ice speed (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(lateral_rate*scyr, 'lateral calving rate (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(calving_dthck, 'calving_dthck (m)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine thickness_based_calving

!---------------------------------------------------------------------------

  subroutine stress_based_calving(&
       nx,                 ny,            &
       dx,                 dy,            &  ! m
       dt,                                &  ! s
       itest,   jtest,     rtest,         &
       calving_front_mask,                &
       speed,                             &  ! m/s
       cf_length,                         &  ! m
       thck_effective,                    &  ! m
       tau_eigen1,    tau_eigen2,         &  ! Pa
       tau_eigenconstant1,                &
       tau_eigenconstant2,                &
       stress_threshold,                  &  ! Pa
       lateral_rate_min,                  &  ! m/s
       calving_dthck)                        ! m

    ! Calve ice based on the eigenvalues of the 2D horizontal stress tensor near the calving front.
    ! When the effective stress is above a certain threshold, the CF advances.
    ! Below the threshold, the CF retreats.

    use glide_diagnostics, only: point_diag

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt                        ! time step (s)

    integer, dimension(nx,ny), intent(in) ::  &
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         speed,                  & ! mean ice speed averaged to cell centers (m/s)
         cf_length,              & ! length of calving front in each grid cell (m)
         thck_effective,         & ! effective thickness for calving (m)
         tau_eigen1, tau_eigen2    ! eigenvalues of the horizontal stress tensor (Pa)

    real(dp), intent(in) :: &
         tau_eigenconstant1,     & ! multiplier for tau_eigen1 (unitless)
         tau_eigenconstant2,     & ! multiplier for tau_eigen2 (unitless)
         stress_threshold,       & ! stress threshold for calving front retreat
         lateral_rate_min          ! min lateral rate at CF (m/s)

    real(dp), dimension(nx,ny), intent(inout) :: &
         calving_dthck             ! thickness lost due to calving (m)

    ! local variables

    integer :: i, j

    real(dp), dimension(nx,ny) :: &
         effec_stress,           & ! effective stress (Pa)
         lateral_rate              ! lateral calving rate (m/s)

    ! Compute thinning in calving-front cells based on the principal eigenvalues of the
    ! horizontal stress tensor, in relation to a stress threshold.
    !
    ! The lateral calving rate Cr is given by
    !
    !   Cr = |v| * effec_stress / stress_threshold
    !   where effec_stress = k1*max(tau_eigen1,0.0) + k2*max(tau_eigen2,0.0)
    !        |v| = ice speed (>=0) at the calving front
    !
    ! In other words:
    ! * Cr = 0 where effec_stress = 0; no calving
    ! * Cr < |v| where effec_stress < stress_threhold
    ! * Cr > |v| where effec_stress > stress_threhold
    ! * Cr = |v| where effec_stress = stress_threhold
    !
    ! The lateral calving rate is converted to a thinning rate dH using
    !
    !   dH = Cr * dt * H_eff * cf_length/ (dx*dy),
    !
    ! The RHS is equal to the ice volume removed per unit grid cell area during one time step

    lateral_rate = 0.0d0
    calving_dthck = 0.0d0
    effec_stress = 0.0d0

    !TODO - Compute nonzero values for floating cells only?
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          effec_stress(i,j) = tau_eigenconstant1 * max(tau_eigen1(i,j), 0.0d0)   &
                            + tau_eigenconstant2 * max(tau_eigen2(i,j), 0.0d0)
       enddo
    enddo

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo

          if (calving_front_mask(i,j) == 1) then
             lateral_rate(i,j)  = speed(i,j) * effec_stress(i,j) / stress_threshold
             lateral_rate(i,j) = max(lateral_rate(i,j), lateral_rate_min)
             calving_dthck(i,j) = lateral_rate(i,j) * dt * thck_effective(i,j) * cf_length(i,j) / (dx*dy)
          endif   ! CF mask

       enddo   ! i
    enddo   ! j

    if (verbose_calving) then
       call point_diag(speed*scyr, 'Stress-based calving, ice speed (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(tau_eigen1, 'tau_eigen1 (Pa)', itest, jtest, rtest, 7, 7, '(f10.0)')
       call point_diag(tau_eigen2, 'tau_eigen2 (Pa)', itest, jtest, rtest, 7, 7, '(f10.0)')
       call point_diag(effec_stress, 'eff_stress (Pa)', itest, jtest, rtest, 7, 7, '(f10.0)')
       call point_diag(lateral_rate*scyr, 'lateral calving rate (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(calving_dthck, 'calving_dthck (m)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine stress_based_calving

!---------------------------------------------------------------------------

  subroutine stochastic_stress_based_calving(&
       nx,                 ny,            &
       dx,                 dy,            &  ! m
       dt,                                &  ! s
       itest,   jtest,     rtest,         &
       parallel,                          &
       calving_front_mask,                &
       thck,                              &  ! m
       thck_effective,                    &  ! m
       effective_areafrac,                &
       speed,                             &  ! m/s
       cf_length,                         &  ! m
       tau_eigen1,    tau_eigen2,         &  ! Pa
       tau_eigenconstant1,                &
       tau_eigenconstant2,                &
       stress_threshold,                  &  ! Pa
       effec_stress_min,                  &  ! Pa
       length_scale,                      &  ! m
       calving_dthck)                        ! m

    ! Calve ice based on the eigenvalues of the 2D horizontal stress tensor near the calving front.
    ! Instead of computing a deterministic calving rate (as in the subroutine above),
    ! let the calving be stochastic.
    !
    ! The general idea is that each cell is assigned a damage probability based on the effective stress.
    ! With high effective stress, it becomes asymptotically more likely to fail.
    ! In each cell with nonzero damage, generate a random number to decide whether it fails.
    ! After identifying failed cells, remove all the failed cells at the CF.
    ! Also remove failed interior cells that have a path through failed cells to the CF.
    !
    ! This idea is based roughly on the percolation model of Bahr (1995).
    ! Still under construction.

    use glissade_masks, only: glissade_fill, initial_color, fill_color, boundary_color
    use glide_diagnostics, only: point_diag

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt                        ! time step (s)

    integer, dimension(nx,ny), intent(in) ::  &
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         thck,                   & ! ice thickness (m)
         thck_effective,         & ! effective ice thickness (m) at the CF
         effective_areafrac,     & !
         speed,                  & ! mean ice speed averaged to cell centers (m/s)
         cf_length,              & ! length of calving front in each grid cell (m)
         tau_eigen1, tau_eigen2    ! eigenvalues of the horizontal stress tensor (Pa)

    real(dp), intent(in) :: &
         tau_eigenconstant1,     & ! multiplier for tau_eigen1 (unitless)
         tau_eigenconstant2,     & ! multiplier for tau_eigen2 (unitless)
         stress_threshold,       & ! stress threshold for calving front retreat (Pa)
         effec_stress_min,       & ! minimum effec stress (Pa) for purposes of computing CF retreat
         length_scale              ! length scale for calving events (m);
                                   ! typically less than or equal to the grid scale

    real(dp), dimension(nx,ny), intent(inout) :: &
         calving_dthck             ! thickness lost due to calving (m)

    ! local variables

    integer :: i, j
    integer :: count, count_save, iter, max_iter

    real(dp), dimension(nx,ny) :: &
         effec_stress,           & ! effective stress (Pa)
         calving_probability       ! probability that all or part of a CF cell is calved during this time step

    real(dp) :: &
         rnd                       ! random number between 0 and 1

    real(dp) :: &
         stress_fac, speed_fac, length_fac  ! unitless factors that determine calving probability

    real(dp) :: grid_scale         ! sqrt(dx*dy) = grid cell length scale (= dx = dy on a square grid)
    real(dp) :: length_ratio       ! calving length_scale / grid_scale

    ! Initialize
    effec_stress = 0.0d0
    calving_probability = 0.0d0
    calving_dthck = 0.0d0
    grid_scale = sqrt(dx*dy)

    ! Compute the effective stress
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          effec_stress(i,j) = tau_eigenconstant1 * max(tau_eigen1(i,j), 0.0d0)   &
                            + tau_eigenconstant2 * max(tau_eigen2(i,j), 0.0d0)
          effec_stress(i,j) = max(effec_stress(i,j), effec_stress_min)
       enddo
    enddo

    if (verbose_calving) then
       call point_diag(tau_eigen1, 'Stochastic calving, tau_eigen1 (Pa)', itest, jtest, rtest, 7, 7, '(f10.0)')
       call point_diag(tau_eigen2, 'tau_eigen2 (Pa)', itest, jtest, rtest, 7, 7, '(f10.0)')
       call point_diag(effec_stress, 'effec_stress (Pa)', itest, jtest, rtest, 7, 7, '(f10.0)')
    endif

    ! Calve cells at the calving front, with a probability based on effec_stress, speed, and other factors.
    ! Notes:
    ! (1) If the timestep is doubled, the probability of a calving event during a given timestep is doubled,
    !      as desired (assuming p << 1). This follows from exp[-(t1 + t2)] = exp(-t1)*exp(-t2).
    ! (2) The calved ice thickness (calving_dthck) is based on thck_effective, not thck.
    !     In CF cells where calving_dthck > thck, the calving will extend into upstream cells.
    !     This is done in subroutine apply_calving_dthck.
    ! (3) If length_scale = grid_scale, then the calving probability p is the probability
    !      of a calving event of extent grid_scale^2 (i.e., one grid cell area) during interval dt.
    !     length_scale < grid_scale, then the calving_probability is greater, but less
    !      area is calved per event. E.g., if length_scale = grid_scale/2, then calving events
    !      are 4x as likely, but the calved area per event is 4x smaller.
    !     In either case, the time-average calving area is the same.
    ! (4) A larger value of length_ratio (i.e., fewer, bigger icebergs) implies greater random variation in CF location.
    ! (5) Deterministic stress-based calving can be viewed as equivalent to stochastic stress-based calving,
    !     in the limit as length_scale -> 0 and there are many small calving events.
    ! (6) The method is probably more realistic for length_ratio < 1 (icebergs smaller than grid scale)
    !     than length_ratio > 1 (icebergs larger than grid scale).

    length_ratio = length_scale/grid_scale

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (calving_front_mask(i,j) == 1) then
             ! compute three unitless factors
             stress_fac = effec_stress(i,j) / stress_threshold
             speed_fac = speed(i,j) * (dt/length_ratio**2) / grid_scale
             length_fac = cf_length(i,j) / grid_scale
             ! compute the calving probability
             calving_probability(i,j) = 1.0d0 - exp(-stress_fac*speed_fac*length_fac)
             call random_number(rnd)
             if (rnd < calving_probability(i,j)) then
                calving_dthck(i,j) = thck_effective(i,j) * length_ratio**2
             endif
          endif   ! CF mask
       enddo   ! i
    enddo   ! j

    call parallel_halo(calving_probability, parallel)
    call parallel_halo(calving_dthck, parallel)

    if (verbose_calving) then
       call point_diag(calving_probability, 'calving probability', itest, jtest, rtest, 7, 7, '(f10.3)')
       call point_diag(calving_dthck, 'calving_dthck (m)', itest, jtest, rtest, 7, 7, '(f10.3)')
    endif

  end subroutine stochastic_stress_based_calving

!---------------------------------------------------------------------------

  subroutine eigencalving(&
       nx,                 ny,            &
       dx,                 dy,            &  ! m
       dt,                                &  ! s
       itest,   jtest,     rtest,         &
       calving_front_mask,                &
       cf_length,                         &  ! m
       thck_effective,                    &  ! m
       eps_eigen1,    eps_eigen2,         &  ! 1/s
       eigenconstant,                     &  ! m
       calving_dthck)                        ! m

    ! Calve ice based on the eigenvalues of the 2D horizontal stress tensor near the calving front,
    !  following Levermann et al. (2012).
    ! Here, eigenconstant corresponds to the Levermann K, but with different units because
    !  it multiplies the geometric mean of strain rates rather than the product.

    use glide_diagnostics, only: point_diag

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt                        ! time step (s)

    integer, dimension(nx,ny), intent(in) ::  &
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         cf_length,              & ! length of calving front in each grid cell (m)
         thck_effective,         & ! effective thickness for calving (m)
         eps_eigen1, eps_eigen2    ! eigenvalues of the horizontal strain rate tensor (1/s)

    real(dp), intent(in) :: &
         eigenconstant             ! strain-rate multiplier (m)

    real(dp), dimension(nx,ny), intent(inout) :: &
         calving_dthck             ! thickness lost due to calving (m)

    ! local variables

    integer :: i, j

    real(dp), dimension(nx,ny) :: &
         lateral_rate               ! lateral calving rate (m/s)

    ! Compute thinning in calving-front cells based on the principal eigenvalues of the
    ! horizontal strain rate tensor.
    !
    ! The lateral calving rate Cr is given by
    !
    !   Cr = K * sqrt[max(0,eps_eigen1 * max(0,eps_eigen2)]
    !
    ! In other words, the calving rate is proportional to the geometric mean of the along-flow
    !  and across-flow strain-rate eigenvalues, provided both eigenvalues are positive.
    !
    ! The lateral calving rate is converted to a thinning rate dH using
    !
    !   dH = Cr * dt * H_eff * cf_length / (dx*dy),
    !
    ! The RHS is equal to the ice volume removed per unit grid cell area during one time step

    lateral_rate = 0.0d0
    calving_dthck = 0.0d0

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (calving_front_mask(i,j) == 1) then

             if (eps_eigen1(i,j) > 0.0d0 .and. eps_eigen2(i,j) > 0.0d0) then
                lateral_rate(i,j) = eigenconstant * sqrt(eps_eigen1(i,j) * eps_eigen2(i,j))
             endif

             calving_dthck(i,j) = lateral_rate(i,j) * dt * thck_effective(i,j) * cf_length(i,j) / (dx*dy)

          endif   ! CF mask
       enddo   ! i
    enddo   ! j

    if (verbose_calving) then
       call point_diag(eps_eigen1*scyr, 'eps_eigen1 (yr^-1)',  itest, jtest, rtest, 7, 7, '(f10.6)')
       call point_diag(eps_eigen2*scyr, 'eps_eigen2 (yr^-1)',  itest, jtest, rtest, 7, 7, '(f10.6)')
       call point_diag(sqrt(max(eps_eigen1,0.d0)*max(eps_eigen2,0.0d0))*scyr, &
            'eps geom mean (yr^-1)',  itest, jtest, rtest, 7, 7, '(f10.6)')
       call point_diag(lateral_rate*scyr, 'lateral calving rate (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(calving_dthck, 'calving_dthck (m)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine eigencalving

!---------------------------------------------------------------------------

  ! This subroutine previously was called for the CALVING_DAMAGE option,
  ! which now uses subroutine stochastic_damage_based_calving.
  ! Keeping the older version for reference and further testing.
  subroutine damage_based_calving(&
       nx,       ny,       nz,            &
       dx,                 dy,            &  ! m
       sigma,              dt,            &
       itest,   jtest,     rtest,         &
       floating_mask,                     &
       calving_front_mask,                &
       speed,                             &  ! m/s
       cf_length,                         &  ! m
       thck_effective,                    &  ! m
       tau_eigen1,    tau_eigen2,         &  ! Pa
       tau_eigenconstant1,                &  ! 1/s
       tau_eigenconstant2,                &  ! 1/s
       damage_threshold,                  &
       damage,                            &
       calving_dthck)                        ! m

    ! Calve ice based on accumulated damage.
    ! This is similar to stress-based calving, except that stresses contribute to damage
    ! instead of directly determining the lateral calving rate.

    use glide_diagnostics, only: point_diag

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny, nz,             & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt                        ! time step (s)

    real(dp), dimension(nz), intent(in) ::&
         sigma                     ! vertical sigma coordinate

    integer, dimension(nx,ny), intent(in) ::  &
         floating_mask,          & ! = 1 where ice is present and floating, else = 0
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         speed,                  & ! mean ice speed averaged to cell centers (m/s)
         cf_length,              & ! length of calving front within a cell (m)
         thck_effective,         & ! effective thickness for calving (m)
         tau_eigen1, tau_eigen2    ! eigenvalues of the horizontal stress tensor (Pa)

    real(dp), intent(in) :: &
         tau_eigenconstant1,     & ! multiplier for tau_eigen1 (unitless)
         tau_eigenconstant2,     & ! multiplier for tau_eigen2 (unitless)
         damage_threshold          ! damage value at which calving begins (unitless)

    real(dp), dimension(nz-1,nx,ny), intent(inout) :: &
         damage                    ! 3D damage tracer, 0 > damage < 1

    real(dp), dimension(nx,ny), intent(inout) :: &
         calving_dthck             ! thickness lost due to calving (m)

    ! local variables

    integer :: i, j, k

    real(dp), dimension(nx,ny) :: &
         damage_column,          & ! average damage in an ice column
         effec_stress,           & ! effective stress (Pa)
         lateral_rate              ! lateral calving rate (m/s)

    real(dp) :: &
         d_damage_dt,            & ! rate of change of damage scalar (1/s)
         damage_frac               ! ratio of damage to damage_threshold (unitless)

    !TODO - Make this a config parameter (similar to damage_threshold, but different units)
    real(dp), parameter :: damage_scale = 1.d7   ! sustained stress over time (Pa yr) required to yield D = 1

    ! Compute thinning in calving-front cells based on damage, a acalar quantity
    !  prognosed from eigenvalues of the horizontal stress tensor (or strain rate tensor).
    !
    ! The lateral calving rate Cr is given by
    !
    !   Cr = |v| *(D/D_thr)
    !   where D = damage in the column
    !         D_thr is a threshold where the calving rate balances the ice speed
    !         |v| = ice speed (>=0) at the calving front
    !
    ! In other words:
    ! * Cr = 0 where D = 0; no calving
    ! * For 0 < D < D_thr, there is some calving, but not enough to halt ice advance.
    ! * For D = D_thr, the calving rate is in balance with the ice speed, so the CF is stable.
    ! * For D > D_thr, the CF retreats.
    !
    ! The lateral calving rate is converted to a thinning rate dH using
    !
    !   dH = Cr * dt * H_eff * cf_length/ (dx*dy),
    !
    ! The RHS is equal to the ice volume removed per unit grid cell area during one time step
    !
    ! Physically, this is similar to the stress-based calving scheme.
    ! The difference is that in the stress-based scheme, the calving rate depends on the stresses at the CF.
    ! In the damage-based scheme, the calving rate depends on the damage at the CF; damage can occur
    !  upstream and then is advected to the CF.

    ! Prognose changes in damage.
    ! For now, this is done using a simple scheme based on the horizontal stress tensor.
    ! At some point, we may want to prognose damage in a way that depends on other factors such as mass balance.
    ! Note: Damage is formally a 3D field, which makes it easier to advect, even though
    !       (in the current scheme) the damage source term is uniform in each column.

    ! Compute the change in damage in each cell
    do j = 2, ny-1
       do i = 1, nx-1
          if (floating_mask(i,j) == 1) then

             effec_stress(i,j) = tau_eigenconstant1 * max(tau_eigen1(i,j), 0.0d0)   &
                               + tau_eigenconstant2 * max(tau_eigen2(i,j), 0.0d0)

             d_damage_dt = effec_stress(i,j) / (damage_scale*scyr)

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
    do j = 2, ny-1
       do i = 2, nx-1
          do k = 1, nz-1
             damage_column(i,j) = damage_column(i,j) + damage(k,i,j) * (sigma(k+1) - sigma(k))
          enddo
       enddo
    enddo

    if (verbose_calving) then
       call point_diag(speed*scyr, 'Damage-based calving, speed (m/yr)', &
            itest, jtest, rtest, 7, 7, '(f10.2)')
       call point_diag(tau_eigen1, 'tau_eigen1 (Pa)',  itest, jtest, rtest, 7, 7, '(f10.0)')
       call point_diag(tau_eigen2, 'tau_eigen2 (Pa)',  itest, jtest, rtest, 7, 7, '(f10.0)')
       call point_diag(effec_stress, 'eff_stress (Pa)',  itest, jtest, rtest, 7, 7, '(f10.0)')
       call point_diag(damage_column, 'new damage', itest, jtest, rtest, 7, 7, '(f10.6)')
    endif

    lateral_rate = 0.0d0
    calving_dthck = 0.0d0

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (calving_front_mask(i,j) == 1) then

             ! Compute a damage fraction, increasing linearly from 0 to 1/damage_threshold.
             damage_frac = damage_column(i,j) / damage_threshold

             ! Convert the damage fraction to a calving rate
             ! For D < D_threshold, the CF will advance.
             ! For D = D_threshold, Cr = speed and the CF should be stable.
             ! For D > D_threshold, the CF will retreat.

             lateral_rate(i,j) = speed(i,j) * damage_frac
             calving_dthck(i,j) = lateral_rate(i,j) * dt * thck_effective(i,j) * cf_length(i,j) / (dx*dy)

          endif   ! CF mask
       enddo   ! i
    enddo   ! j

    if (verbose_calving) then
       call point_diag(lateral_rate*scyr, 'lateral calving rate (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(calving_dthck, 'calving_dthck (m)', itest, jtest, rtest, 7, 7)
    endif

  end subroutine damage_based_calving

!---------------------------------------------------------------------------

  subroutine stochastic_damage_based_calving(&
       nx,       ny,       nz,            &
       dx,                 dy,            &  ! m
       sigma,              dt,            &
       itest,   jtest,     rtest,         &
       parallel,                          &
       floating_mask,                     &
       calving_front_mask,                &
       thck,                              &  ! m
       topg,                              &  ! m
       tau_eigen1,    tau_eigen2,         &  ! Pa
       tau_eigenconstant1,                &  ! 1/s
       tau_eigenconstant2,                &  ! 1/s
       stress_threshold,                  &  ! Pa
       effec_stress_min,                  &  ! Pa
       damage_constant,                   &  ! Pa s
       damage,                            &
       calving_dthck,                     &  ! m
       eps_eigen1,    eps_eigen2)            ! 1/s

    ! Calve ice based on accumulated damage.
    ! This is similar to stress-based calving, except that stresses contribute to damage
    ! instead of directly determining the lateral calving rate.

    use glide_diagnostics, only: point_diag
    use glissade_masks, only: glissade_fill, initial_color, fill_color, boundary_color

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny, nz,             & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt                        ! time step (s)

    real(dp), dimension(nz), intent(in) ::&
         sigma                     ! vertical sigma coordinate

    integer, dimension(nx,ny), intent(in) ::  &
         floating_mask,          & ! = 1 where ice is present and floating, else = 0
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         thck,                   & ! ice thickness (m)
         topg,                   & ! bed topography (m)
         tau_eigen1, tau_eigen2    ! eigenvalues of the horizontal stress tensor (Pa)

    real(dp), intent(in) :: &
         tau_eigenconstant1,     & ! multiplier for tau_eigen1 (unitless)
         tau_eigenconstant2,     & ! multiplier for tau_eigen2 (unitless)
         stress_threshold,       & ! stress value at which damage begins (Pa)
         effec_stress_min,       & ! min value for effec_stress (Pa)
         damage_constant           ! damage scaling term (Pa s)

    real(dp), dimension(nz-1,nx,ny), intent(inout) :: &
         damage                    ! 3D damage tracer, 0 > damage < 1

    real(dp), dimension(nx,ny), intent(inout) :: &
         calving_dthck             ! thickness lost due to calving (m)

    real(dp), dimension(nx,ny), intent(in), optional :: &
         eps_eigen1, eps_eigen2    ! eigenvalues of the horizontal strain rate tensor (s^-1)

    ! local variables

    integer :: i, j, k

    integer, dimension(nx,ny) :: &
         damage_mask               ! = 1 for damaged cells, else = 0

    real(dp), dimension(nx,ny) :: &
         damage_column,          & ! average damage in an ice column
         effec_stress              ! effective stress (Pa)

    real(dp) :: &
         d_damage_dt,            & ! rate of change of damage scalar (1/s)
         rnd                       ! random number

    integer :: count, count_save, iter, max_iter

    integer,  dimension(nx,ny) ::  &
         color                     ! integer 'color' for identifying cells that will calve

    logical, parameter :: tau2_only = .false.

    ! Compute calving based on damage, a scalar quantity prognosed from eigenvalues
    ! of the horizontal stress tensor (or strain rate tensor).
    !
    ! This is a stochastic version of a damage-based scheme. The idea is as follows:
    ! * Damage is a scalar D in the range [0,1]. At each timestep, the pre-existing damage field is advected with the ice.
    ! * In this subroutine, first compute changes in damage due to stress.
    ! * The resulting field D can have any values in the range [0,1].
    ! * Use a random number to convert damage into a binary field, with values of 0 and 1 only.
    ! * Use a flood-fill algorithm to remove connected regions of damaged ice (D = 1) at the CF.

    ! Compute the change in damage in each cell
    ! For now, this is done using a simple scheme based on the horizontal stress tensor.
    ! At some point, we may want to prognose damage in a way that depends on other factors such as mass balance.
    ! Note: Damage is formally a 3D field, which makes it easier to advect, even though
    !       (in the current scheme) the damage source term is uniform in each column.

    ! Compute an effective stress for damage, based on the eigenvalues of the horizontal stress tensor

    effec_stress = 0.0d0

    do j = 2, ny-1
       do i = 1, nx-1
          if (floating_mask(i,j) == 1) then
             effec_stress(i,j) = &
                  max(tau_eigenconstant1*tau_eigen1(i,j) + tau_eigenconstant2*tau_eigen2(i,j), 0.0d0)
          endif
       enddo
    enddo

    ! In calving-front cells, impose a lower bound for effec_stress
    do j = 2, ny-1
       do i = 1, nx-1
          if (calving_front_mask(i,j) == 1) then
             effec_stress(i,j) = max(effec_stress(i,j), effec_stress_min)
          endif
       enddo
    enddo

    ! Compute the change in damage as a function of effec_stress
    do j = 2, ny-1
       do i = 1, nx-1
          if (floating_mask(i,j) == 1) then
             !WHL - Assume a threshold of zero?
!!             d_damage_dt = (effec_stress(i,j) - stress_threshold) / damage_constant
             ! Note: damage_constant has units of Pa*s, so d_damage_dt has units of 1/s
             d_damage_dt = effec_stress(i,j) / damage_constant
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
    do j = 2, ny-1
       do i = 2, nx-1
          do k = 1, nz-1
             damage_column(i,j) = damage_column(i,j) + damage(k,i,j) * (sigma(k+1) - sigma(k))
          enddo
       enddo
    enddo

    call parallel_halo(damage_column, parallel)

    if (verbose_calving) then
       call point_diag(tau_eigen1, 'tau_eigen1 (Pa)',  itest, jtest, rtest, 9, 9, '(f8.0)')
       call point_diag(tau_eigen2, 'tau_eigen2 (Pa)',  itest, jtest, rtest, 9, 9, '(f8.0)')
       call point_diag(effec_stress, 'effec_stress (Pa)',  itest, jtest, rtest, 9, 9, '(f8.0)')
       call point_diag(sqrt(tau_eigen1**2 + tau_eigen2**2), 'true eff stress (Pa)',  itest, jtest, rtest, 9, 9, '(f8.0)')
       call point_diag(damage_column, 'new damage', itest, jtest, rtest, 9, 9, '(f8.4)')
    endif

    ! Convert damage to a binary field, with D = 0.0 or 1.0 everywhere

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (damage_column(i,j) > 0.0d0) then
             call random_number(rnd)
             if (rnd < damage_column(i,j)) then
                damage_column(i,j) = 1.0d0
                damage(:,i,j) = 1.0d0
             else
                damage_column(i,j) = 0.0d0
                damage(:,i,j) = 0.0d0
             endif
          endif   ! damage > 0
       enddo   ! i
    enddo   ! j

    call parallel_halo(damage_column, parallel)

    ! Based on damage_column, compute a binary integer mask
    where (damage_column == 1.0d0)
       damage_mask = 1
    elsewhere
       damage_mask = 0
    endwhere

    if (verbose_calving) then
       call point_diag(damage_mask, 'damage_mask', itest, jtest, rtest, 9, 9, '(i8)')
    endif

    ! Mark all floating cells with the initial color.
    ! Mark other cells with the boundary color.

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (floating_mask(i,j) == 1) then
             color(i,j) = initial_color
          else
             color(i,j) = boundary_color
          endif
       enddo
    enddo

    call parallel_halo(color, parallel)

    ! Seed the fill with damaged CF cells, and then recursively fill all damaged cells that
    !  have a path to the CF through other damaged cells and are therefore to be calved.
    ! Repeat the recursion as necessary to spread the fill to adjacent processors.

    max_iter = max(parallel%ewtasks, parallel%nstasks)
    count_save = 0

    do iter = 1, max_iter

       if (iter == 1) then

          do j = 1, ny
             do i = 1, nx
                if (calving_front_mask(i,j) == 1 .and. damage_mask(i,j) == 1) then  ! cell that can seed the fill
                   if (color(i,j) /= fill_color) then
                      ! assign the fill color to this cell, and recursively fill neighbor cells
                      call glissade_fill(&
                           nx,    ny,    &
                           i,     j,     &
                           color, damage_mask)
                   endif
                endif
             enddo
          enddo

       else   ! iter > 1

          call parallel_halo(color, parallel)

          ! Loop through halo cells that were just filled on neighbor processors

          ! west halo layer
          i = nhalo
          do j = 1, ny
             if (color(i,j) == fill_color) then
                call glissade_fill(&
                     nx,    ny,    &
                     i+1,   j,     &
                     color, damage_mask)
             endif
          enddo

          ! east halo layer
          i = nx - nhalo + 1
          do j = 1, ny
             if (color(i,j) == fill_color) then
                call glissade_fill(&
                     nx,    ny,    &
                     i-1,   j,     &
                     color, damage_mask)
             endif
          enddo

          ! south halo layer
          j = nhalo
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color) then
                call glissade_fill(&
                     nx,    ny,    &
                     i,     j+1,   &
                     color, damage_mask)
             endif
          enddo

          ! north halo layer
          j = ny-nhalo+1
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color) then
                call glissade_fill(&
                     nx,    ny,    &
                     i,     j-1,   &
                     color, damage_mask)
             endif
          enddo

       endif   ! iter

       ! Count the calving cells
       count = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (color(i,j) == fill_color) count = count + 1
          enddo
       enddo
       !WHL - If running a large problem, may want to reduce the frequency of this global sum
       count = parallel_reduce_sum(count)

       if (count == count_save) then
!!          if (verbose_calving .and. main_task) &
!!               write(6,*) 'Fill converged: iter, global_count =', iter, global_count
          exit
       else
!!          if (verbose_calving .and. main_task) &
!!               write(6,*) 'Convergence check: iter, global_count =', iter, global_count
          count_save = count
       endif

    enddo   ! max_iter

    call parallel_halo(color, parallel)

    ! Once the fill is done, any cells with the fill color are to be calved.
    ! The actual calving is done in another subroutine.

    where (color == fill_color .and. floating_mask == 1)
       calving_dthck = thck
    elsewhere
       calving_dthck = 0.0d0
    endwhere

    if (verbose_calving) then
!       if (main_task) then
!          print*, 'boundary_color, initial color, fill color:', &
!               boundary_color, initial_color, fill_color
!       endif
       call point_diag(color, 'After damage, color', itest, jtest, rtest, 9, 9, '(i8)')
       call point_diag(calving_dthck, 'calving_dthck (m)', itest, jtest, rtest, 9, 9, '(f8.0)')
    endif

  end subroutine stochastic_damage_based_calving

!---------------------------------------------------------------------------

  subroutine calving_front_advance_retreat(&
       nx,                 ny,             &
       dx,                 dy,             &
       dt,                 time,           &  ! s
       itest,   jtest,     rtest,          &
       calving_front_mask,                 &
       thck_pre_transport,                 &  ! m
       thck,                               &  ! m
       cf_length,                          &  ! m
       thck_effective,                     &  ! m
       cf_advance_retreat_amplitude,       &  ! m/s
       cf_advance_retreat_period,          &  ! s
       calving_dthck)                         ! m

    ! Calve ice based on a prescribed rate of calving front advance or retreat.
    ! This option requires a subgrid calving front scheme.

    use glimmer_physcon, only: pi
    use glide_diagnostics, only: point_diag

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt,                     & ! time step (s)
         time                      ! elapsed time (s) of model run

    integer, dimension(nx,ny), intent(in)  ::  &
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         thck_pre_transport        ! ice thickness (m) before doing transport, SMB, and BMB;
                                   ! determines the calving needed to keep the CF stationary

    real(dp), intent(in) :: &
         cf_advance_retreat_amplitude, & ! prescribed amplitude (m/yr) for calving front advance or retreat
                                         ! positive for sin(2*pi*t/period), negative for -sin(2*pi*t/period)
                                         ! should be *negative* for CalvingMIP Experiments 2 and 4
         cf_advance_retreat_period       ! period (s) for an advance/retreat cycle
                                         ! period = 0 => constant amplitude


    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m)
         cf_length,              & ! length of calving front in each grid cell (m)
         thck_effective            ! effective thickness for calving (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         calving_dthck             ! thickness (m) to be added to the calving flux

    ! local variables

    integer :: i, j

    real(dp), dimension(nx,ny) :: &
         dthck_dt_transport        ! rate of thickness change from transport, SMB and BMB

    real(dp) :: &
         cf_advance_retreat_rate,& ! prescribed rate of CF advance or retreat (m/s); positive for advance
         cf_net_advance            ! total advance (or if negative, retreat) since time 0

    ! Set the calving front advance/retreat rate
    if (cf_advance_retreat_period > 0.0d0) then
       cf_advance_retreat_rate = &
            cf_advance_retreat_amplitude * sin(2.d0*pi*time / cf_advance_retreat_period)
       cf_net_advance = cf_advance_retreat_amplitude * (cf_advance_retreat_period/(2.d0*pi)) &
            * (1.d0 - cos(2.d0*pi*time/cf_advance_retreat_period))
    else
       cf_advance_retreat_rate = cf_advance_retreat_amplitude ! > 0 for advance, < 0 for retreat
    endif

    dthck_dt_transport = (thck - thck_pre_transport) / dt
    calving_dthck = 0.0d0

    if (verbose_calving) then
       call point_diag(dthck_dt_transport*scyr, 'dH/dt transport (m/yr)', itest, jtest, rtest, 7, 7)
    endif

    if (verbose_calving .and. this_rank==rtest) then
       write(6,*) ' '
       write(6,*) 'Time (yr)', time/scyr
       write(6,*) 'CF advance_retreat_rate (m/yr) =', cf_advance_retreat_rate*scyr
       write(6,*) 'Net CF advance (km):', cf_net_advance / 1000.d0
       write(6,*) 'Prescribed radius (km):', 750.d0 + cf_net_advance/1000.d0
       write(6,*) 'Prescribed circumference (km):', 2.d0*pi*(750.d0 + cf_net_advance/1000.d0)
    endif

    ! Loop over locally owned cells
    ! Calving occurs only in CF cells: ice-filled cells with one or more ocean neighbors.

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo

          if (calving_front_mask(i,j) == 1) then

             ! First, compute the calving needed to compensate for the thickness increase in partial CF cells.
             ! This step should hold the calving front stationary by balancing transport/SMB/BMB with calving.
             ! Note the sign: calving_dthck > 0 for a positive calving flux, with thinning in cell (i,j).

             calving_dthck(i,j) = max(dthck_dt_transport(i,j), 0.0d0) * dt

             ! If the CF is supposed to retreat, then increase the calving (negative increment for dthck).
             ! If the CF is supposed to advance, then reduce the calving (positive increment for dthck).
             ! Note: Some calving might already have been done in an unprotected cell just past the current CF.
             !       If so, it is possible to undo this calving so the CF can advance.

             ! Decrease the calving if the CF is advancing, increase if the CF is retreating

             calving_dthck(i,j) = calving_dthck(i,j) - &
                  (cf_advance_retreat_rate*dt * thck_effective(i,j) * cf_length(i,j)) / (dx*dy)

             ! Make sure dthck >= 0. This prevents unphysical CF advance via negative calving
             calving_dthck(i,j) = max(calving_dthck(i,j), 0.0d0)

          endif   ! calving_front cell
       enddo   ! i
    enddo   ! j

  end subroutine calving_front_advance_retreat

!HG<---------------------------------------------------------------------------

  subroutine calving_front_meltrate(&
       nx,                 ny,             &
       dx,                 dy,             &
       dt,                 time,           &  ! s
       itest,   jtest,     rtest,          &
       calving_front_mask,                 &
       thck_pre_transport,                 &  ! m
       thck,                               &  ! m
       cf_length,                          &  ! m
       thck_effective,                     &  ! m
       cf_advance_retreat_amplitude,       &  ! m/s
       cf_advance_retreat_period,          &  ! s
       calving_dthck)                         ! m

    ! Melt ice "horizontally" based on a prescribed melt rate of equvivalent calving front retreat.
    ! This option requires a subgrid calving front scheme.

    use glimmer_physcon, only: pi
    use glide_diagnostics, only: point_diag

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt,                     & ! time step (s)
         time                      ! elapsed time (s) of model run

    integer, dimension(nx,ny), intent(in)  ::  &
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         thck_pre_transport        ! ice thickness (m) before doing transport, SMB, and BMB;
                                   ! determines the calving needed to keep the CF stationary

    real(dp), intent(in) :: &
         cf_advance_retreat_amplitude, & ! prescribed amplitude (m/s) for calving front advance or retreat
                                         ! positive for sin(2*pi*t/period), negative for -sin(2*pi*t/period)
                                         ! should be *negative* for CalvingMIP Experiments 2 and 4
         cf_advance_retreat_period       ! period (s) for an advance/retreat cycle
                                         ! period = 0 => constant amplitude


    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m)
         cf_length,              & ! length of calving front in each grid cell (m)
         thck_effective            ! effective thickness for calving (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         calving_dthck             ! thickness (m) to be added to the calving flux

    ! local variables

    integer :: i, j

    real(dp), dimension(nx,ny) :: &
         dthck_dt_transport        ! rate of thickness change from transport, SMB and BMB

    real(dp) :: &
         cf_advance_retreat_rate,& ! prescribed rate of CF advance or retreat (m/s); positive for advance
         cf_net_advance            ! total advance (or if negative, retreat) since time 0

    ! Set the calving front advance/retreat rate
    if (cf_advance_retreat_period > 0.0d0) then
       cf_advance_retreat_rate = &
            cf_advance_retreat_amplitude * sin(2.d0*pi*time / cf_advance_retreat_period)
       cf_net_advance = cf_advance_retreat_amplitude * (cf_advance_retreat_period/(2.d0*pi)) &
            * (1.d0 - cos(2.d0*pi*time/cf_advance_retreat_period))
    else
       cf_advance_retreat_rate = cf_advance_retreat_amplitude ! > 0 for advance, < 0 for retreat
    endif

    dthck_dt_transport = (thck - thck_pre_transport) / dt
    calving_dthck = 0.0d0

    if (verbose_calving) then
       call point_diag(dthck_dt_transport*scyr, 'dH/dt transport (m/yr)', itest, jtest, rtest, 7, 7)
    endif

    if (verbose_calving .and. this_rank==rtest) then
       write(6,*) ' '
       write(6,*) 'Time (yr)', time/scyr
       write(6,*) 'CF meltrate (m/yr) =', cf_advance_retreat_rate*scyr
       write(6,*) 'Net CF advance (km):', cf_net_advance / 1000.d0
       write(6,*) 'Prescribed radius (km):', 750.d0 + cf_net_advance/1000.d0
       write(6,*) 'Prescribed circumference (km):', 2.d0*pi*(750.d0 + cf_net_advance/1000.d0)
    endif

    ! Loop over locally owned cells
    ! Calving occurs only in CF cells: ice-filled cells with one or more ocean neighbors.

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo

          if (calving_front_mask(i,j) == 1) then

!HG< commented to apply an absolute “horizontal” melt rate - i.e., the melt rate perpendicular to the face of an approximately vertical calving front             
!             ! First, compute the calving needed to compensate for the thickness increase in partial CF cells.
!             ! This step should hold the calving front stationary by balancing transport/SMB/BMB with calving.
!             ! Note the sign: calving_dthck > 0 for a positive calving flux, with thinning in cell (i,j).
!
!             calving_dthck(i,j) = max(dthck_dt_transport(i,j), 0.0d0) * dt

             ! If the CF is supposed to retreat, then increase the calving (negative increment for dthck).
             ! If the CF is supposed to advance, then reduce the calving (positive increment for dthck).
             ! Note: Some calving might already have been done in an unprotected cell just past the current CF.
             !       If so, it is possible to undo this calving so the CF can advance.

             ! Decrease the calving if the CF is advancing, increase if the CF is retreating

             calving_dthck(i,j) = calving_dthck(i,j) - &
                  (cf_advance_retreat_rate*dt * thck_effective(i,j) * cf_length(i,j)) / (dx*dy)

             ! Make sure dthck >= 0. This prevents unphysical CF advance via negative calving
             calving_dthck(i,j) = max(calving_dthck(i,j), 0.0d0)
             
             if (verbose_calving) then
                write(6,*) ' '
                write(6,*) 'Time (yr)', time/scyr
                write(6,*) 'applied melt as calving_dthck (m) =', i,j, calving_dthck(i,j), thck(i,j), thck_effective(i,j), cf_length(i,j)
             endif

          endif   ! calving_front cell
       enddo   ! i
    enddo   ! j

  end subroutine calving_front_meltrate

!HG<---------------------------------------------------------------------------

  subroutine slater_calving_meltrate(&
       nx,                 ny,             &
       dx,                 dy,             &
       dt,                 time,           &  ! s
       itest,   jtest,     rtest,          &
       calving_front_mask,                 &
       thck_pre_transport,                 &  ! m
       thermal_forcing_applied,            &  ! degC
       runoff_applied,                     &  ! m/s
       thck,                               &  ! m
       cf_length,                          &  ! m
       thck_effective,                     &  ! m
       cf_advance_retreat_amplitude,       &  ! m/s
       cf_advance_retreat_period,          &  ! s
       calving_dthck)                         ! m

    ! This routine does two things 
    ! 1. Calving based on Slater/Rahlves c = alpha * Heff^beta
    ! 2. Melt ice "horizontally" based on a prescribed melt rate of equvivalent calving front retreat.
    ! Melt parameterised as a function of runoff and thermal_forcing, both 2d fields
    ! This option requires a subgrid calving front scheme.

    use glimmer_physcon, only: pi
    use glide_diagnostics, only: point_diag

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt,                     & ! time step (s)
         time                      ! elapsed time (s) of model run

    integer, dimension(nx,ny), intent(in)  ::  &
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         thck_pre_transport        ! ice thickness (m) before doing transport, SMB, and BMB;
                                   ! determines the calving needed to keep the CF stationary
    real(dp), dimension(nx,ny), intent(in) :: &
         thermal_forcing_applied,& ! thermal_forcing applied (degC)
         runoff_applied            ! runoff_applied (m/s)
    
    real(dp), intent(in) :: &
         cf_advance_retreat_amplitude, & ! prescribed amplitude (m/s) for calving front advance or retreat
                                         ! positive for sin(2*pi*t/period), negative for -sin(2*pi*t/period)
                                         ! should be *negative* for CalvingMIP Experiments 2 and 4
         cf_advance_retreat_period       ! period (s) for an advance/retreat cycle
                                         ! period = 0 => constant amplitude


    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m)
         cf_length,              & ! length of calving front in each grid cell (m)
         thck_effective            ! effective thickness for calving (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         calving_dthck             ! thickness (m) to be added to the calving flux

    ! local variables

    integer :: i, j

    real(dp), dimension(nx,ny) :: &
         dthck_dt_transport        ! rate of thickness change from transport, SMB and BMB

    real(dp) :: &
         alpha_sr, beta_sr,&       ! parameters in Slater/Rahlves calving law
         c_sr,&                    ! calving rate in m/yr calculated from Slater/Rahlves law
         q_sr,&                    ! runoff in Rignot calculation in m/d.  
         tf_sr,&                   ! thermal forcing in deg C
         m_sr,&                    ! melting rate in m/yr calculated from Slater ISMIP6 melt approach
         cf_advance_retreat_rate,& ! prescribed rate of CF advance or retreat (m/s); positive for advance
         cf_net_advance            ! total advance (or if negative, retreat) since time 0

    ! Set the calving front advance/retreat rate
    if (cf_advance_retreat_period > 0.0d0) then
       cf_advance_retreat_rate = &
            cf_advance_retreat_amplitude * sin(2.d0*pi*time / cf_advance_retreat_period)
       cf_net_advance = cf_advance_retreat_amplitude * (cf_advance_retreat_period/(2.d0*pi)) &
            * (1.d0 - cos(2.d0*pi*time/cf_advance_retreat_period))
    else
       cf_advance_retreat_rate = cf_advance_retreat_amplitude ! > 0 for advance, < 0 for retreat
    endif

    dthck_dt_transport = (thck - thck_pre_transport) / dt
    calving_dthck = 0.0d0

    if (verbose_calving) then
       call point_diag(dthck_dt_transport*scyr, 'dH/dt transport (m/yr)', itest, jtest, rtest, 7, 7)
    endif

    if (verbose_calving .and. this_rank==rtest) then
       write(6,*) ' '
       write(6,*) 'Time (yr)', time/scyr
       write(6,*) 'CF meltrate (m/yr) =', cf_advance_retreat_rate*scyr
       write(6,*) 'Net CF advance (km):', cf_net_advance / 1000.d0
       write(6,*) 'Prescribed radius (km):', 750.d0 + cf_net_advance/1000.d0
       write(6,*) 'Prescribed circumference (km):', 2.d0*pi*(750.d0 + cf_net_advance/1000.d0)
    endif

    ! Loop over locally owned cells
    ! Calving occurs only in CF cells: ice-filled cells with one or more ocean neighbors.

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo

          if (calving_front_mask(i,j) == 1) then

             ! 1. Apply Slater/Rahlves calving law c = alpha * Heff^beta, assuming beta=3 and alpha=1.5e-5
             ! For a 700 m calving front, this should give a calving front retreat rate of ~ 5000 m/yr

             alpha_sr = 0.5e-5
             beta_sr = 3.0
             c_sr = -alpha_sr * thck_effective(i,j)**beta_sr/scyr ! m/s
             
             ! Increase the 'calving' if the CF is retreating/melting
             calving_dthck(i,j) = calving_dthck(i,j) - &
                  (c_sr*dt * thck_effective(i,j) * cf_length(i,j)) / (dx*dy)

             if (verbose_calving .and. thck(i,j) > 0.1) then
                print*, 'slater calving'
                print*, i, j, thck_effective(i,j), thck(i,j), alpha_sr, beta_sr, c_sr, calving_dthck(i,j)
             endif

             ! 2. Apply an absolute “horizontal” melt rate - i.e., the melt rate perpendicular to the face of an approximately vertical calving front             

             ! Implemented as a global melt rate using negative values for config parameter cf_advance_retreat_amplitude
             ! Increase the 'calving' if the CF is retreating/melting

             !calving_dthck(i,j) = calving_dthck(i,j) - &
             !     (cf_advance_retreat_rate*dt * thck_effective(i,j) * cf_length(i,j)) / (dx*dy)

             ! Or parametersing melt as function of runoff, thermal forcing and effective thickness
             ! See Rignot et al. 2016 or ISMIP6 melt forcing approach.

             !tf_sr = 2.0 ! hard coding a cold cavity thermal forcing of 2 degC
             !q_sr = 10.0 ! hard coding a runoff of 400 m3/s multiplied by s/day and divided by cross-sectional area of 700m*5000m [m/day] 
             
             tf_sr = thermal_forcing_applied(i,j) ! 2d thermal forcing [degC]
             q_sr = runoff_applied(i,j) * 86400.  ! runoff_applied passed in m/s; for Rignot equation convert to [m/d] 

             m_sr = -(3.0 * 1.0e-4 * thck_effective(i,j) * q_sr**0.39 + 0.15) * tf_sr**1.18 * 365./scyr ! Rignot et al. 2016; formulted in m/d, converted to m/s
             
             calving_dthck(i,j) = calving_dthck(i,j) - &
                  (m_sr*dt * thck_effective(i,j) * cf_length(i,j)) / (dx*dy)

             if (verbose_calving .and. thck(i,j) > 0.1) then
                print*, 'slater melting'
                print*, i, j, thck_effective(i,j), thck(i,j), tf_sr, q_sr, m_sr, calving_dthck(i,j)
             endif
             ! Make sure dthck >= 0. This prevents unphysical CF advance via negative calving
             calving_dthck(i,j) = max(calving_dthck(i,j), 0.0d0)

          endif   ! calving_front cell
       enddo   ! i
    enddo   ! j

  end subroutine slater_calving_meltrate

!HG<---------------------------------------------------------------------------

  subroutine slater_meltrate(&
       nx,                 ny,             &
       dx,                 dy,             &
       dt,                 time,           &  ! s
       itest,   jtest,     rtest,          &
       calving_front_mask,                 &
       thck_pre_transport,                 &  ! m
       thermal_forcing_applied,            &  ! degC
       runoff_applied,                     &  ! m/s
       thck,                               &  ! m
       cf_length,                          &  ! m
       thck_effective,                     &  ! m
       cf_advance_retreat_amplitude,       &  ! m/s
       cf_advance_retreat_period,          &  ! s
       calving_dthck)                         ! m

    ! Melt ice "horizontally" based on a prescribed melt rate of equvivalent calving front retreat.
    ! Melt parameterised as a function of runoff and thermal_forcing, both 2d fields
    ! This option requires a subgrid calving front scheme.

    use glimmer_physcon, only: pi
    use glide_diagnostics, only: point_diag

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy,                 & ! grid cell size (m)
         dt,                     & ! time step (s)
         time                      ! elapsed time (s) of model run

    integer, dimension(nx,ny), intent(in)  ::  &
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(nx,ny), intent(in) :: &
         thck_pre_transport        ! ice thickness (m) before doing transport, SMB, and BMB;
                                   ! determines the calving needed to keep the CF stationary
    real(dp), dimension(nx,ny), intent(in) :: &
         thermal_forcing_applied,& ! thermal_forcing applied (degC)
         runoff_applied            ! runoff_applied (m/s)
    
    real(dp), intent(in) :: &
         cf_advance_retreat_amplitude, & ! prescribed amplitude (m/s) for calving front advance or retreat
                                         ! positive for sin(2*pi*t/period), negative for -sin(2*pi*t/period)
                                         ! should be *negative* for CalvingMIP Experiments 2 and 4
         cf_advance_retreat_period       ! period (s) for an advance/retreat cycle
                                         ! period = 0 => constant amplitude


    real(dp), dimension(nx,ny), intent(inout) :: &
         thck,                   & ! ice thickness (m)
         cf_length,              & ! length of calving front in each grid cell (m)
         thck_effective            ! effective thickness for calving (m)

    real(dp), dimension(nx,ny), intent(out) :: &
         calving_dthck             ! thickness (m) to be added to the calving flux

    ! local variables

    integer :: i, j

    real(dp), dimension(nx,ny) :: &
         dthck_dt_transport        ! rate of thickness change from transport, SMB and BMB

    real(dp) :: &
         q_sr,&                    ! runoff in Rignot calculation in m/d.  
         tf_sr,&                   ! thermal forcing in deg C
         m_sr,&                    ! melting rate in m/yr calculated from Slater ISMIP6 melt approach
         cf_advance_retreat_rate,& ! prescribed rate of CF advance or retreat (m/s); positive for advance
         cf_net_advance            ! total advance (or if negative, retreat) since time 0

    ! Set the calving front advance/retreat rate
    if (cf_advance_retreat_period > 0.0d0) then
       cf_advance_retreat_rate = &
            cf_advance_retreat_amplitude * sin(2.d0*pi*time / cf_advance_retreat_period)
       cf_net_advance = cf_advance_retreat_amplitude * (cf_advance_retreat_period/(2.d0*pi)) &
            * (1.d0 - cos(2.d0*pi*time/cf_advance_retreat_period))
    else
       cf_advance_retreat_rate = cf_advance_retreat_amplitude ! > 0 for advance, < 0 for retreat
    endif

    dthck_dt_transport = (thck - thck_pre_transport) / dt
    calving_dthck = 0.0d0

    if (verbose_calving) then
       call point_diag(dthck_dt_transport*scyr, 'dH/dt transport (m/yr)', itest, jtest, rtest, 7, 7)
    endif

    if (verbose_calving .and. this_rank==rtest) then
       write(6,*) ' '
       write(6,*) 'Time (yr)', time/scyr
       write(6,*) 'CF meltrate (m/yr) =', cf_advance_retreat_rate*scyr
       write(6,*) 'Net CF advance (km):', cf_net_advance / 1000.d0
       write(6,*) 'Prescribed radius (km):', 750.d0 + cf_net_advance/1000.d0
       write(6,*) 'Prescribed circumference (km):', 2.d0*pi*(750.d0 + cf_net_advance/1000.d0)
    endif

    ! Loop over locally owned cells
    ! Calving occurs only in CF cells: ice-filled cells with one or more ocean neighbors.

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo

          if (calving_front_mask(i,j) == 1) then

             ! Apply an absolute “horizontal” melt rate - i.e., the melt rate perpendicular to the face of an approximately vertical calving front             

             ! Could be implemented as a global melt rate using negative values for config parameter cf_advance_retreat_amplitude
             ! Increase the 'calving' if the CF is retreating/melting

             !calving_dthck(i,j) = calving_dthck(i,j) - &
             !     (cf_advance_retreat_rate*dt * thck_effective(i,j) * cf_length(i,j)) / (dx*dy)

             ! Here parametersing melt as function of runoff, thermal forcing and effective thickness
             ! See Rignot et al. 2016 or ISMIP6 melt forcing approach.

             !tf_sr = 2.0 ! hard coding a cold cavity thermal forcing of 2 degC
             !q_sr = 10.0 ! hard coding a runoff of 400 m3/s multiplied by s/day and divided by cross-sectional area of 700m*5000m [m/day] 
             
             tf_sr = thermal_forcing_applied(i,j) ! 2d thermal forcing [degC]
             q_sr = runoff_applied(i,j) * 86400.  ! runoff_applied passed in m/s; for Rignot equation convert to [m/d] 

             m_sr = -(3.0 * 1.0e-4 * thck_effective(i,j) * q_sr**0.39 + 0.15) * tf_sr**1.18 * 365./scyr ! Rignot et al. 2016; formulted in m/d, converted to m/s
             
             calving_dthck(i,j) = calving_dthck(i,j) - &
                  (m_sr*dt * thck_effective(i,j) * cf_length(i,j)) / (dx*dy)

             if (verbose_calving .and. thck(i,j) > 0.1) then
                print*, 'slater melting'
                print*, i, j, thck_effective(i,j), thck(i,j), tf_sr, q_sr, m_sr, calving_dthck(i,j)
             endif
             ! Make sure dthck >= 0. This prevents unphysical CF advance via negative calving
             calving_dthck(i,j) = max(calving_dthck(i,j), 0.0d0)

          endif   ! calving_front cell
       enddo   ! i
    enddo   ! j

  end subroutine slater_meltrate

!---------------------------------------------------------------------------

  subroutine apply_calving_dthck(&
       nx,           ny,       &
       itest, jtest, rtest,    &
       parallel,               &
       calving_front_mask,     &
       floating_mask,          &
       full_mask,              &
       flux_in,                &
       calving_dthck,          &
       thck,                   &
       calving_thck)

    ! Apply calving_dthck as computed from a given calving law.
    ! This is the thinning rate that will give the desired lateral calving rate.

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    integer, dimension(nx,ny), intent(in)  ::  &
         calving_front_mask,     & ! = 1 where ice is floating and borders at least one ocean cell, else = 0
         floating_mask,          & ! = 1 where ice is present and floating, else = 0
         full_mask                 ! = 1 for ice-filled cells that are not partial_cf cells, else = 0

    real(dp), dimension(-1:1,-1:1,nx,ny), intent(in) :: &
         flux_in                   ! ice volume fluxes (m^3/s) into cell from each neighbor cell

    real(dp), dimension(nx,ny), intent(inout) :: &
         calving_dthck,          & ! thickness (m) to be added to the calving flux
         thck,                   & ! ice thickness (m)
         calving_thck              ! thickness lost due to calving (m)

    ! local variables

    integer :: i, j, ii, jj, iup, jup
    integer :: count
    real(dp) :: total_flux, total_dthck, my_dthck

    ! Apply calving_dthck to calving-front cells

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (calving_front_mask(i,j) == 1) then

             my_dthck = calving_dthck(i,j)
             if (my_dthck > thck(i,j)) then
                ! calve the full column and compute the remainder available for upstream calving
                calving_dthck(i,j) = calving_dthck(i,j) - thck(i,j)
                calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                thck(i,j) = 0.0d0
             else
                ! calve part of the column
                thck(i,j) = thck(i,j) - my_dthck
                calving_thck(i,j) = calving_thck(i,j) + my_dthck
                calving_dthck(i,j) = 0.0d0
             endif

          endif   ! CF mask
       enddo   ! i
    enddo   ! j

    call parallel_halo(thck, parallel)
    call parallel_halo(calving_dthck, parallel)

    ! If calving_dthck > 0 after a full column is calved, apply further calving to upstream neighbors.
    ! The thinning in each neighbor is proportional to that cell's contribution to the downstream flux.
    ! Loop includes halo cells so each local cell can be thinned by all its downstream neighbors.

    do j = 2, ny-1
       do i = 2, nx-1
          if (calving_dthck(i,j) > 0.0d0) then

             ! Given flux_in (ice flux in m^3/s entering the cell from each upstream neighbor),
             ! compute the fraction of the thinning to be applied to each upstream neighbor.
             count = 0
             total_flux = 0.0d0
             do jj = -1, 1
                do ii = -1, 1
                   iup = i + ii; jup = j + jj
                   !TODO - Is this logic correct? Both floating and full?
                   if (flux_in(ii,jj,i,j) > 0.0d0 .and. floating_mask(iup,jup) == 1 .and. &
                        full_mask(iup,jup) == 1) then
                      count = count + 1
                      total_flux = total_flux + flux_in(ii,jj,i,j)
                   endif
                enddo
             enddo

             !WHL - debug
             if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                print*, 'Calve upstream: dthck, input flux (m^3/yr)=', calving_dthck(i,j), total_flux*scyr
                print*, '   No. of upstream cells =', count
             endif

             ! Calve ice in the upstream neighbors
             if (total_flux > 0.0d0) then
                total_dthck = calving_dthck(i,j)
                do jj = -1, 1
                   do ii = -1, 1
                      iup = i + ii; jup = j + jj
                      if (flux_in(ii,jj,i,j) > 0.0d0 .and. floating_mask(iup,jup) == 1 &
                           .and. full_mask(iup,jup) == 1) then
                         my_dthck = total_dthck * flux_in(ii,jj,i,j)/total_flux
                         thck(iup,jup) = thck(iup,jup) - my_dthck
                         calving_thck(iup,jup) = calving_thck(iup,jup) + my_dthck
                         calving_dthck(i,j) = calving_dthck(i,j) - my_dthck
!                         if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
!                            print*, '   Upstream ii, jj, frac, remaining calving_dthck:', &
!                                 ii, jj, flux_in(ii,jj,i,j)/total_flux, calving_dthck(i,j)
!                         endif
                      endif
                   enddo
                enddo
             endif

          endif   ! calving_dthck > 0
       enddo   ! i
    enddo   ! j

  end subroutine apply_calving_dthck

!---------------------------------------------------------------------------

  subroutine apply_grounded_calving_dthck(&
       nx,           ny,       &
       itest, jtest, rtest,    &
       parallel,               &
       calving_front_mask,     &
       ice_mask,               &
       floating_mask,          &
       full_mask,              &
       flux_in,                &
       calving_dthck,          &
       thck,                   &
       calving_thck)

    ! Apply calving_dthck as computed from a given calving law.
    ! This is the thinning rate that will give the desired lateral calving rate.

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    integer, dimension(nx,ny), intent(in)  ::  &
         calving_front_mask,     & ! = 1 where ice is floating and borders at least one ocean cell, else = 0
         ice_mask,               & ! = 1 where ice is present, else = 0
         floating_mask,          & ! = 1 where ice is present and floating, else = 0
         full_mask                 ! = 1 for ice-filled cells that are not partial_cf cells, else = 0

    real(dp), dimension(-1:1,-1:1,nx,ny), intent(in) :: &
         flux_in                   ! ice volume fluxes (m^3/s) into cell from each neighbor cell

    real(dp), dimension(nx,ny), intent(inout) :: &
         calving_dthck,          & ! thickness (m) to be added to the calving flux
         thck,                   & ! ice thickness (m)
         calving_thck              ! thickness lost due to calving (m)

    ! local variables

    integer :: i, j, ii, jj, iup, jup
    integer :: count
    real(dp) :: total_flux, total_dthck, my_dthck

    ! Apply calving_dthck to calving-front cells

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (calving_front_mask(i,j) == 1) then

             my_dthck = calving_dthck(i,j)
             if (my_dthck > thck(i,j)) then
                ! calve the full column and compute the remainder available for upstream calving
                calving_dthck(i,j) = calving_dthck(i,j) - thck(i,j)
                calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                thck(i,j) = 0.0d0
             else
                ! calve part of the column
                thck(i,j) = thck(i,j) - my_dthck
                calving_thck(i,j) = calving_thck(i,j) + my_dthck
                calving_dthck(i,j) = 0.0d0
             endif

          endif   ! CF mask
       enddo   ! i
    enddo   ! j

    call parallel_halo(thck, parallel)
    call parallel_halo(calving_dthck, parallel)

    ! If calving_dthck > 0 after a full column is calved, apply further calving to upstream neighbors.
    ! The thinning in each neighbor is proportional to that cell's contribution to the downstream flux.
    ! Loop includes halo cells so each local cell can be thinned by all its downstream neighbors.

    do j = 2, ny-1
       do i = 2, nx-1
          if (calving_dthck(i,j) > 0.0d0) then

             ! Given flux_in (ice flux in m^3/s entering the cell from each upstream neighbor),
             ! compute the fraction of the thinning to be applied to each upstream neighbor.
             count = 0
             total_flux = 0.0d0
             do jj = -1, 1
                do ii = -1, 1
                   iup = i + ii; jup = j + jj
                   !TODO - Is this logic correct? Both grounded and full?
                   if (flux_in(ii,jj,i,j) > 0.0d0 .and. ice_mask(iup,jup) == 1 &
                        .and. floating_mask(iup,jup) == 0 .and. full_mask(iup,jup) == 1) then
                      count = count + 1
                      total_flux = total_flux + flux_in(ii,jj,i,j)
                   endif
                enddo
             enddo

             !WHL - debug
             if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                print*, 'Calve upstream: dthck, input flux (m^3/yr)=', calving_dthck(i,j), total_flux*scyr
                print*, '   No. of upstream cells =', count
             endif

             ! Calve ice in the upstream neighbors
             if (total_flux > 0.0d0) then
                total_dthck = calving_dthck(i,j)
                do jj = -1, 1
                   do ii = -1, 1
                      iup = i + ii; jup = j + jj
                      !TODO - Is this logic correct? Both grounded and full?
                      if (flux_in(ii,jj,i,j) > 0.0d0 .and. ice_mask(iup,jup) == 1 &
                           .and. floating_mask(iup,jup) == 0 .and. full_mask(iup,jup) == 1) then
                         my_dthck = total_dthck * flux_in(ii,jj,i,j)/total_flux
                         thck(iup,jup) = thck(iup,jup) - my_dthck
                         calving_thck(iup,jup) = calving_thck(iup,jup) + my_dthck
                         calving_dthck(i,j) = calving_dthck(i,j) - my_dthck
                         if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                            print*, '   Upstream ii, jj, frac, remaining calving_dthck:', &
                                 ii, jj, flux_in(ii,jj,i,j)/total_flux, calving_dthck(i,j)
                         endif
                      endif
                   enddo
                enddo
             endif

          endif   ! calving_dthck > 0
       enddo   ! i
    enddo   ! j

  end subroutine apply_grounded_calving_dthck

!---------------------------------------------------------------------------

  subroutine advance_calving_front(&
       nx,           ny,     &
       itest, jtest, rtest,  &
       ocean_mask,           &
       calving_front_mask,   &
       flux_in,              &
       thck_effective,       &
       thck)

    ! Check for thck > thck_effective in CF cells. This can happen if ice in unprotected cells
    ! has been redistributed to full or nearly full cells upstream, and then was not calved.
    ! Distribute excess ice downstream.

    ! input/output arguments

    integer, intent(in) :: &
         nx, ny                    ! grid dimensions

    integer, intent(in) :: &
         itest, jtest, rtest       ! coordinates of diagnostic point

    integer, dimension(nx,ny), intent(in)  ::  &
         ocean_mask,             & ! = 1 where topg is below sea level and ice is absent, else = 0
         calving_front_mask        ! = 1 where ice is floating and borders at least one ocean cell, else = 0

    real(dp), dimension(-1:1,-1:1,nx,ny), intent(in) :: &
         flux_in                   ! ice volume fluxes (m^3/s) into cell from each neighbor cell

    real(dp), dimension(nx,ny), intent(in) :: &
         thck_effective            ! effective thickness (m) at the CF

    real(dp), dimension(nx,ny), intent(inout) :: &
         thck                      ! ice thickness (m)

    ! local variables

    integer :: i, j, ii, jj, idn, jdn
    integer :: count
    real(dp) :: total_flux, total_dthck, my_dthck

    ! Omitting this call is not answer-changing
!!    call parallel_halo(ocean_mask, parallel)

    ! Loop must include halo cells so that each local cell can receive contributions
    !  from all its upstream neighbors
    do j = 2, ny-1
       do i = 2, nx-1
          if (thck(i,j) > thck_effective(i,j) .and. calving_front_mask(i,j) == 1) then

             ! Compute the fraction of excess ice in this cell to give to each downstream ocean neighbor
             count = 0
             total_flux = 0.0d0
             do jj = -1, 1
                do ii = -1, 1
                   idn = i + ii; jdn = j + jj
                   !TODO - Is this logic OK?
                   if (flux_in(-ii,-jj,idn,jdn) > 0.0d0 .and. ocean_mask(idn,jdn) == 1) then
                      count = count + 1
                      total_flux = total_flux + flux_in(-ii,-jj,idn,jdn)
                   endif
                enddo
             enddo

             if (verbose_calving .and. abs(i-itest)<=1 .and. abs(j-jtest)<=1 .and. this_rank==rtest) then
                print*, ' '
                print*, 'Excess ice: rank, i, j, dthck, downstream flux (m^3/s)=', &
                     this_rank, i, j, thck(i,j) - thck_effective(i,j), total_flux
                print*, '   No. of downstream cells =', count
             endif

             ! Move ice to its downstream ocean neighbors
             total_dthck = thck(i,j) - thck_effective(i,j)
             if (total_flux > 0.0d0) then
                do jj = -1, 1
                   do ii = -1, 1
                      idn = i + ii; jdn = j + jj
                      if (flux_in(-ii,-jj,idn,jdn) > 0.0d0 .and. ocean_mask(idn,jdn) == 1) then
                         my_dthck = total_dthck * flux_in(-ii,-jj,idn,jdn)/total_flux
                         thck(idn,jdn) = thck(idn,jdn) + my_dthck
                         thck(i,j) = thck(i,j) - my_dthck
                         if (verbose_calving .and. abs(i-itest)<=1 .and. abs(j-jtest)<=1 .and. this_rank==rtest) then
                            print*, '   Downstream ii, jj, frac:', ii, jj, flux_in(-ii,-jj,idn,jdn)/total_flux
                         endif
                      endif
                   enddo
                enddo
             endif

          endif   ! CF cell with thck > thck_effective
       enddo   ! i
    enddo   ! j

  end subroutine advance_calving_front

!---------------------------------------------------------------------------
  
  subroutine glissade_remove_icebergs(&
       nx,           ny,            &
       parallel,                    &
       itest, jtest, rtest,         &
       f_ground_threshold,          &
       thck,                        &
       f_ground_cell,               &
       ice_mask,                    &
       floating_mask,               &
       land_mask,                   &
       calving_thck)

    ! Remove any icebergs
        
    ! The algorithm is as follows:
    ! (1) Mark all cells with ice with the initial color.
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
    !Note: When using a subgrid CF scheme, it is safer to pass in full_mask in place of the usual ice_mask,
    !      so that partial CF cells do not spread the fill.
    integer,  dimension(nx,ny), intent(inout) :: ice_mask        !> = 1 where ice is present (thck > thklim), else = 0;
                                                                 !> may exclude partial CF cells
    integer,  dimension(nx,ny), intent(in)    :: floating_mask   !> = 1 where ice is present and floating, else = 0
    integer,  dimension(nx,ny), intent(in)    :: land_mask       !> = 1 where topg - eus >= 0, else = 0
    real(dp), dimension(nx,ny), intent(inout) :: calving_thck    !> thickness lost due to calving in each grid cell;
                                                                 !> on output, includes ice in icebergs
    ! local variables

    integer :: i, j, iter

    integer :: &
         max_iter,             & ! max(ewtasks, nstasks)
         local_count,          & ! local counter for filled values
         global_count,         & ! global counter for filled values
         global_count_save       ! global counter for filled values from previous iteration

    integer,  dimension(nx,ny) ::  &
         color                 ! integer 'color' for identifying icebergs

    if (verbose_calving) then
       call point_diag(thck, 'Remove icebergs, thck (m)', itest, jtest, rtest, 7, 7)
       call point_diag(ice_mask, 'ice_mask', itest, jtest, rtest, 7, 7)
!!       call point_diag(f_ground_cell, 'f_ground_cell', itest, jtest, rtest, 7, 7)
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

                if (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0 .and. &
                    f_ground_cell(i,j) >= f_ground_threshold) then  ! grounded ice

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

       else  ! iter > 1

          ! Check for halo cells that were just filled on neighbor processors
          ! TODO - Is the 'ice_mask = 1' logic redundant?
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

          ! east halo layer
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

       endif  ! iter

       local_count = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (color(i,j) == fill_color) local_count = local_count + 1
          enddo
       enddo

       !WHL - If running a large problem, may want to reduce the frequency of this global sum
       global_count = parallel_reduce_sum(local_count)

       if (global_count == global_count_save) then
!!          if (verbose_calving .and. main_task) &
!!               write(6,*) 'Fill converged: iter, global_count =', iter, global_count
          exit
       else
!!          if (verbose_calving .and. main_task) &
!!               write(6,*) 'Convergence check: iter, global_count =', iter, global_count
          global_count_save = global_count
       endif

    enddo  ! max_iter

    ! Icebergs are cells that still have the initial color and are not on land.
    ! Remove ice in these cells, adding it to the calving field.

    do j = 2, ny-1
       do i = 2, nx-1
          if (color(i,j) == initial_color .and. land_mask(i,j) == 0) then
             calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
             thck(i,j) = 0.0d0
             !TODO - Also handle tracers?  E.g., set damage(:,i,j) = 0.d0?
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
    !  or thin floating ice (H < 10 m) on two opposite sides (i.e., in cells (i-1,j) and (i+1,j),
    !  or (i,j-1) and (i,j+1)), connected to ice on the other two sides.
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

    real(dp) :: mask_e, mask_w, mask_n, mask_s
    real(dp) :: thin_threshold

    ! Both floating and weakly grounded cells can be identified as isthmuses and removed;
    !  f_ground_threshold is used to identify weakly grounded cells.

    ! An isthmus cell has ice-free ocean or thin floating ice on each side:
    !  ice is considered thin and floating if thinner than isthmus_thck_threshold
    !  and also thinner than the ice in the central cell.

    real(dp), parameter :: &   ! threshold (m) for counting floating ice as thin
         isthmus_thck_threshold = 10.0d0

    if (verbose_calving) then
       call point_diag(thck, 'Remove isthmuses, thck', itest, jtest, rtest, 7, 7)
    endif

    do j = 2, ny-1
       do i = 2, nx-1
          if (floating_mask(i,j) == 1 .or. f_ground_cell(i,j) < f_ground_threshold) then
             mask_e = 0; mask_w = 0; mask_n = 0; mask_s = 0
             thin_threshold = min(isthmus_thck_threshold, thck(i,j))
             if (ocean_mask(i+1,j) == 1 .or. &
                  (floating_mask (i+1,j) == 1 .and. thck(i+1,j) < thin_threshold) ) then
                mask_e = 1
             endif
             if (ocean_mask(i-1,j) == 1 .or. &
                  (floating_mask (i-1,j) == 1 .and. thck(i-1,j) < thin_threshold) ) then
                mask_w = 1
             endif
             if (ocean_mask(i,j+1) == 1 .or. &
                  (floating_mask (i,j+1) == 1 .and. thck(i,j+1) < thin_threshold) ) then
                mask_n = 1
             endif
             if (ocean_mask(i,j-1) == 1 .or. &
                  (floating_mask (i,j-1) == 1 .and. thck(i,j-1) < thin_threshold) ) then
                mask_s = 1
             endif
             ! Note: We call a cell an isthmus only if it borders floating ice on the other two edges.
             !       Without this requirement, we will remove many cells that are not true isthmuses.
             if (mask_e == 1 .and. mask_w == 1) then
                if (floating_mask(i,j-1) == 1 .and. floating_mask (i,j+1) == 1) then
                   calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                   thck(i,j) = 0.0d0
                endif
             elseif (mask_n == 1 .and. mask_s == 1) then
                if (floating_mask(i-1,j) == 1 .and. floating_mask (i+1,j) == 1) then
                   calving_thck(i,j) = calving_thck(i,j) + thck(i,j)
                   thck(i,j) = 0.0d0
                endif
             endif
          endif  ! floating or lightly grounded
       enddo  ! i
    enddo  ! j

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

  subroutine extrapolate_to_calving_front(&
       nx,              ny,             &
       partial_cf_mask, full_mask,      &
       effective_areafrac,              &
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

    real(dp), dimension(nx,ny), intent(in) :: &
         effective_areafrac

    real(dp), dimension(nx,ny), intent(inout) :: &
         field1                ! field with values to be extrapolated to partial CF cells

    real(dp), dimension(nx,ny), intent(inout), optional :: &
         field2                ! second field with values to be extrapolated to partial CF cells

    ! local variables

    integer :: i, j
    real(dp) :: field_max

    ! For each partial CF cell, set the field to a weighted combination of the local value
    ! and the value from an adjacent full cell. The local value gets a weight of effective_areafrac.
    ! TODO: Try taking the average instead of the max.

    do j = 2, ny-1
       do i = 2, nx-1
          if (partial_cf_mask(i,j) == 1) then
             field_max = max(field1(i-1,j)*full_mask(i-1,j), field1(i+1,j)*full_mask(i+1,j), &
                             field1(i,j-1)*full_mask(i,j-1), field1(i,j+1)*full_mask(i,j+1))
             field1(i,j) = effective_areafrac(i,j)  * field1(i,j) +  &
                  (1.0d0 - effective_areafrac(i,j)) * field_max
             if (present(field2)) then
                field_max = max(field2(i-1,j)*full_mask(i-1,j), field2(i+1,j)*full_mask(i+1,j), &
                                field2(i,j-1)*full_mask(i,j-1), field2(i,j+1)*full_mask(i,j+1))
                field2(i,j) = effective_areafrac(i,j)  * field2(i,j) +  &
                     (1.0d0 - effective_areafrac(i,j)) * field_max
             endif
          endif
       enddo
    enddo

  end subroutine extrapolate_to_calving_front

!---------------------------------------------------------------------------
! The next two subroutines are diagnostic subroutines for CalvingMIP.
! They estimate the calving front location along 8 prescribed axes
!  for the circular and Thule domains.
! They are not necessary if we have offline tools for locating the CF,
!  but are left here for reference.
!---------------------------------------------------------------------------

  subroutine locate_calving_front_circular(&
       nx,             ny,           &
       dx,             dy,           &
       x0,             y0,           &
       x1,             y1,           &
       parallel,                     &
       itest, jtest,   rtest,        &
       areafrac,                     &
       cf_location)

    use cism_parallel, only: parallel_reduce_maxloc, parallel_reduce_minloc, broadcast
    use glissade_grid_operators, only: glissade_stagger
    use glide_diagnostics, only: point_diag

    ! Find the calving front location along eight profiles on the circular domain.
    ! These profiles are the four cardinal directions (N, S, E, W) along with the diagonals
    !  that form 45-degree angles with the cardinal directions.

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy                    ! grid cell size (m)

    real(dp), dimension(nx-1), intent(in) :: x0  ! x coordinate of NE cell corners
    real(dp), dimension(ny-1), intent(in) :: y0  ! y coordinate of NE cell corners
    real(dp), dimension(nx), intent(in) :: x1  ! x coordinate of cell centers
    real(dp), dimension(ny), intent(in) :: y1  ! y coordinate of cell centers

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    real(dp), dimension(nx,ny), intent(in) :: areafrac
    real(dp), dimension(2,8), intent(out) :: cf_location

    ! local variables

    integer :: i, j, iglobal, jglobal
    integer :: axis
    integer :: procnum
    real(dp) :: cf_location_xmax, cf_location_ymax, cf_location_xmin, cf_location_ymin, radius
    real(dp) :: &
         this_areafrac_avg, next_areafrac_avg  ! average of areafrac in two adjacent cells
    real(dp) :: areafrac_ne, areafrac_nw, areafrac_se, areafrac_sw

    ! Note: For the original CalvingMIP grid, the origin was located at a cell center,
    !        so both axes passed through cell centers.
    !       For the new CalvingMIP grid (as of Nov. 2024), the origin is located at a cell corner,
    !        so both axes lie along cell edges.
    !       The following code is written generally to find the CF along the x-axis and y-axis
    !        in either case.
    !       The code aborts if one of the other is not true.

    logical :: &
         x_axis_thru_centers,  & ! true if the x-axis passes through cell centers
         x_axis_thru_edges,    & ! true if the x-axis passes through cell edges
         y_axis_thru_centers,  & ! true if the y-axis passes through cell centers
         y_axis_thru_edges       ! true if the y-axis passes through cell edges

    ! Find the x and y coordinates of the calving front along the different axes
    !  specified in CalvingMIP.
    ! The code assumes a circular domain with center at (0,0).
    ! The logic depends on whether the N, S, E and W axes pass through cell centers or edges.

    if (this_rank == rtest) print*, 'Locate_calving_front for calvingMIP, rtest =', rtest

    ! Determine whether the x and y axes passes through cell centers, or through cell edges.
    ! They should pass through one or the other.
    x_axis_thru_centers = .false.
    do j = nhalo+1, ny-nhalo
       if (y1(j) == 0.0d0) then
          x_axis_thru_centers = .true.
       endif
    enddo

    x_axis_thru_edges = .false.
    do j = nhalo+1, ny-nhalo
       if (y0(j) == 0.0d0) then
          x_axis_thru_edges = .true.
       endif
    enddo

    y_axis_thru_centers = .false.
    do i = nhalo+1, nx-nhalo
       if (x1(i) == 0.0d0) then
          y_axis_thru_centers = .true.
       endif
    enddo

    y_axis_thru_edges = .false.
    do i = nhalo+1, nx-nhalo
       if (x0(i) == 0.0d0) then
          y_axis_thru_edges = .true.
       endif
    enddo

!    if (x_axis_thru_centers) then
!       print*, this_rank, 'x_axis_thru_centers', x_axis_thru_centers
!    endif
!    if (y_axis_thru_centers) then
!       print*, this_rank, 'y_axis_thru_centers', y_axis_thru_centers
!    endif
!    if (x_axis_thru_edges) then
!       print*, this_rank, 'x_axis_thru_edges', x_axis_thru_edges
!    endif
!    if (y_axis_thru_edges) then
!       print*, this_rank, 'y_axis_thru_edges', y_axis_thru_edges
!    endif

    cf_location(:,:) = 0.0d0

    ! Find the CF location along each of 8 axes
    ! The CF lies in the last cell along a given axis with areafrac > 0
    ! All loops are over locally owned cells

    axis = 1  ! index for the positive y-axis (profile A)
    if (y_axis_thru_centers) then
       do i = nhalo+1, nx-nhalo
          if (x1(i) == 0.0d0) then
             cf_location(1,axis) = 0.0d0
             do j = nhalo+1, ny-nhalo
                if (areafrac(i,j) > 0.0d0 .and. areafrac(i,j+1) == 0.0d0) then
                   cf_location(2,axis) = y1(j) + (areafrac(i,j) - 0.5d0)*dy
                endif
             enddo
          endif
       enddo
    elseif (y_axis_thru_edges) then
       do i = nhalo+1, nx-nhalo
          if (x0(i) == 0.0d0) then  ! E edge of cell lies on the y-axis
             cf_location(1,axis) = 0.0d0
             do j = nhalo+1, ny-nhalo
                this_areafrac_avg = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
                next_areafrac_avg = 0.5d0 * (areafrac(i,j+1) + areafrac(i+1,j+1))
                if (this_areafrac_avg > 0.0d0 .and. next_areafrac_avg == 0.0d0) then
                   cf_location(2,axis) = y1(j) + (this_areafrac_avg - 0.5d0)*dy
                endif
             enddo
          endif
       enddo
    endif   ! y_axis_thru_centers

    ! If this proc has a positive value of y, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_location(2,axis), xout=cf_location_ymax, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 2  ! index for the line y = x in the positive x and y direction (profile B)
    do i = nhalo+1, nx-nhalo
       do j = nhalo+1, ny-nhalo
          if (x1(i) == y1(j)) then ! on the line y = x
             if (areafrac(i,j) > 0.0d0 .and. areafrac(i+1,j+1) == 0.0d0) then
                areafrac_ne = 0.5d0 * (areafrac(i+1,j) + areafrac(i,j+1))
                areafrac_sw = 0.5d0 * (areafrac(i,j-1) + areafrac(i-1,j))
                if (areafrac_ne >= 0.5d0) then  ! CF in cell (i+1,j+1)
                   cf_location(1,axis) = x0(i) + (areafrac_ne - 0.5d0)/areafrac_ne * (0.5d0*dx)
                   cf_location(2,axis) = y0(j) + (areafrac_ne - 0.5d0)/areafrac_ne * (0.5d0*dy)
                elseif (areafrac_sw < 0.5d0) then  ! CF in cell (i-1,j-1)
                   cf_location(1,axis) = x0(i-1) - (0.5d0 - areafrac_sw)/(1.0d0 - areafrac_sw) * (0.5d0*dx)
                   cf_location(2,axis) = y0(j-1) - (0.5d0 - areafrac_sw)/(1.0d0 - areafrac_sw) * (0.5d0*dy)
                else   ! CF in cell (i,j)
                   if (areafrac(i,j) >= 0.5d0) then   ! CF in upper right of cell
                      cf_location(1,axis) = x1(i) + (areafrac(i,j) - 0.5d0)/(areafrac(i,j) - areafrac_ne) * (0.5d0*dx)
                      cf_location(2,axis) = y1(j) + (areafrac(i,j) - 0.5d0)/(areafrac(i,j) - areafrac_ne) * (0.5d0*dy)
                   else  ! areafrac(i,j) < 0.5; CF in lower left of cell
                      cf_location(1,axis) = x1(i) - (0.5d0 - areafrac(i,j))/(areafrac_sw - areafrac(i,j)) * (0.5d0*dx)
                      cf_location(2,axis) = y1(j) - (0.5d0 - areafrac(i,j))/(areafrac_sw - areafrac(i,j)) * (0.5d0*dy)
                   endif
                endif
             endif
          endif   ! on the line y = x
       enddo   ! i
    enddo   ! j

    ! If this proc has a positive value of x, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_location(1,axis), xout=cf_location_ymax, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 3  ! index for the positive x-axis (profile C)
    if (x_axis_thru_centers) then
       do j = nhalo+1, ny-nhalo
          if (y1(j) == 0.0d0) then
             cf_location(2,axis) = 0.0d0
             do i = nhalo+1, nx-nhalo
                if (areafrac(i,j) > 0.0d0 .and. areafrac(i+1,j) == 0.0d0) then
                   cf_location(1,axis) = x1(i) + (areafrac(i,j) - 0.5d0)*dx
                endif
             enddo
          endif
       enddo
    elseif (x_axis_thru_edges) then
       do j = nhalo+1, ny-nhalo
          if (y0(j) == 0.0d0) then
             cf_location(2,axis) = 0.0d0
             do i = nhalo+1, nx-nhalo
                this_areafrac_avg = 0.5d0 * (areafrac(i,j) + areafrac(i,j+1))
                next_areafrac_avg = 0.5d0 * (areafrac(i+1,j) + areafrac(i+1,j+1))
                if (this_areafrac_avg > 0.0d0 .and. next_areafrac_avg == 0.0d0) then
                   cf_location(1,axis) = x1(i) + (this_areafrac_avg - 0.5d0)*dx
                endif
             enddo
          endif
       enddo
    endif   ! x_axis_thru_centers

    ! If this proc has a positive value of x, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_location(1,axis), xout=cf_location_xmax, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 4  ! index for the line y = -x in the positive x and negative y direction (profile D)
    do i = nhalo+1, nx-nhalo
       do j = nhalo+1, ny-nhalo
          if (x1(i) == -y1(j)) then ! on the line y = -x
             if (areafrac(i,j) > 0.0d0 .and. areafrac(i+1,j-1) == 0.0d0) then
                areafrac_se = 0.5d0 * (areafrac(i+1,j) + areafrac(i,j-1))
                areafrac_nw = 0.5d0 * (areafrac(i-1,j) + areafrac(i,j+1))
                if (areafrac_se >= 0.5d0) then  ! CF in cell (i+1,j-1)
                   cf_location(1,axis) = x0(i) + (areafrac_se - 0.5d0)/areafrac_se * (0.5d0*dx)
                   cf_location(2,axis) = y0(j-1) - (areafrac_se - 0.5d0)/areafrac_se * (0.5d0*dy)
                elseif (areafrac_nw < 0.5d0) then  ! CF in cell (i-1,j+1)
                   cf_location(1,axis) = x0(i-1) - (0.5d0 - areafrac_nw)/(1.0d0 - areafrac_nw) * (0.5d0*dx)
                   cf_location(2,axis) = y0(j) + (0.5d0 - areafrac_nw)/(1.0d0 - areafrac_nw) * (0.5d0*dy)
                else   ! CF in cell (i,j)
                   if (areafrac(i,j) >= 0.5d0) then   ! CF in lower right of cell
                      cf_location(1,axis) = x1(i) + (areafrac(i,j) - 0.5d0)/(areafrac(i,j) - areafrac_se) * (0.5d0*dx)
                      cf_location(2,axis) = y1(j) - (areafrac(i,j) - 0.5d0)/(areafrac(i,j) - areafrac_se) * (0.5d0*dy)
                   else  ! areafrac(i,j) < 0.5; CF in upper left of cell
                      cf_location(1,axis) = x1(i) - (0.5d0 - areafrac(i,j))/(areafrac_nw - areafrac(i,j)) * (0.5d0*dx)
                      cf_location(2,axis) = y1(j) + (0.5d0 - areafrac(i,j))/(areafrac_nw - areafrac(i,j)) * (0.5d0*dy)
                   endif
                endif
             endif
          endif
       enddo   ! i
    enddo   ! j

    ! If this proc has a positive value of x, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_location(1,axis), xout=cf_location_ymax, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 5  ! index for the negative y-axis (profile E)
    if (y_axis_thru_centers) then
       do i = nhalo+1, nx-nhalo
          if (x1(i) == 0.0d0) then
             cf_location(1,axis) = 0.0d0
             do j = ny-nhalo, nhalo+1, -1
                if (areafrac(i,j) > 0.0d0 .and. areafrac(i,j-1) == 0.0d0) then
                   cf_location(2,axis) = y1(j) + (0.5d0 - areafrac(i,j))*dy
                endif
             enddo
          endif
       enddo
    elseif (y_axis_thru_edges) then
       do i = nhalo+1, nx-nhalo
          if (x0(i) == 0.0d0) then  ! E edge of cell lies on the y-axis
             cf_location(1,axis) = 0.0d0
             do j = ny-nhalo, nhalo+1, -1
                this_areafrac_avg = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
                next_areafrac_avg = 0.5d0 * (areafrac(i,j-1) + areafrac(i+1,j-1))
                if (this_areafrac_avg > 0.0d0 .and. next_areafrac_avg == 0.0d0) then
                   cf_location(2,axis) = y1(j) + (0.5d0 - this_areafrac_avg)*dy
                endif
             enddo
          endif
       enddo
    endif   ! y_axis_thru_centers

    ! If this proc has a negative value of y, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_location(2,axis), xout=cf_location_ymin, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 6  ! index for the line y = x in the negative x and y direction (profile F)
    do i = nhalo+1, nx-nhalo
       do j = nhalo+1, ny-nhalo
          if (x1(i) == y1(j)) then ! on the line y = x
             if (areafrac(i,j) > 0.0d0 .and. areafrac(i-1,j-1) == 0.0d0) then
                areafrac_sw = 0.5d0 * (areafrac(i,j-1) + areafrac(i-1,j))
                areafrac_ne = 0.5d0 * (areafrac(i,j+1) + areafrac(i+1,j))
                if (areafrac_sw >= 0.5d0) then  ! CF in cell (i-1,j-1)
                   cf_location(1,axis) = x0(i-1) - (areafrac_sw - 0.5d0)/areafrac_sw * (0.5d0*dx)
                   cf_location(2,axis) = y0(j-1) - (areafrac_sw - 0.5d0)/areafrac_sw * (0.5d0*dy)
                elseif (areafrac_ne < 0.5d0) then  ! CF in cell (i+1,j+1)
                   cf_location(1,axis) = x0(i) + (0.5d0 - areafrac_ne)/(1.0d0 - areafrac_ne) * (0.5d0*dx)
                   cf_location(2,axis) = y0(j) + (0.5d0 - areafrac_ne)/(1.0d0 - areafrac_ne) * (0.5d0*dy)
                else   ! CF in cell (i,j)
                   if (areafrac(i,j) >= 0.5d0) then   ! CF in lower left of cell
                      cf_location(1,axis) = x1(i) - (areafrac(i,j) - 0.5d0)/(areafrac(i,j) - areafrac_sw) * (0.5d0*dx)
                      cf_location(2,axis) = y1(j) - (areafrac(i,j) - 0.5d0)/(areafrac(i,j) - areafrac_sw) * (0.5d0*dy)
                   else  ! areafrac(i,j) < 0.5; CF in upper right of cell
                      cf_location(1,axis) = x1(i) + (0.5d0 - areafrac(i,j))/(areafrac_ne - areafrac(i,j)) * (0.5d0*dx)
                      cf_location(2,axis) = y1(j) + (0.5d0 - areafrac(i,j))/(areafrac_ne - areafrac(i,j)) * (0.5d0*dy)
                   endif
                endif
             endif
          endif
       enddo   ! i
    enddo   ! j

    ! If this proc has a negative value of x, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_location(1,axis), xout=cf_location_ymax, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 7  ! index for the negative x-axis (profile g)
    if (x_axis_thru_centers) then
       do j = nhalo+1, ny-nhalo
          if (y1(j) == 0.0d0) then
             cf_location(2,axis) = 0.0d0
             do i = nx-nhalo, nhalo+1, -1
                if (areafrac(i,j) > 0.0d0 .and. areafrac(i-1,j) == 0.0d0) then
                   cf_location(1,axis) = x1(i) + (0.5d0 - areafrac(i,j))*dx
                endif
             enddo
          endif
       enddo
    elseif (x_axis_thru_edges) then
       do j = nhalo+1, ny-nhalo
          if (y0(j) == 0.0d0) then
             cf_location(2,axis) = 0.0d0
             do i = nx-nhalo, nhalo+1, -1
                this_areafrac_avg = 0.5d0 * (areafrac(i,j) + areafrac(i,j+1))
                next_areafrac_avg = 0.5d0 * (areafrac(i-1,j) + areafrac(i-1,j+1))
                if (this_areafrac_avg > 0.0d0 .and. next_areafrac_avg == 0.0d0) then
                   cf_location(1,axis) = x1(i) + (0.5d0 - this_areafrac_avg)*dx
                endif
             enddo
          endif
       enddo
    endif   ! x_axis_thru_centers

    ! If this proc has a negative value of x, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_location(1,axis), xout=cf_location_xmin, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 8  ! index for the line y = -x in the negative x and positive y direction (profile H)
    do i = nhalo+1, nx-nhalo
       do j = nhalo+1, ny-nhalo
          if (x1(i) == -y1(j)) then ! on the line y = -x
             if (areafrac(i,j) > 0.0d0 .and. areafrac(i-1,j+1) == 0.0d0) then
                areafrac_nw = 0.5d0 * (areafrac(i-1,j) + areafrac(i,j+1))
                areafrac_se = 0.5d0 * (areafrac(i+1,j) + areafrac(i,j-1))
                if (areafrac_nw >= 0.5d0) then  ! CF in cell (i-1,j+1)
                   cf_location(1,axis) = x0(i-1) - (areafrac_nw - 0.5d0)/areafrac_nw * (0.5d0*dx)
                   cf_location(2,axis) = y0(j)   + (areafrac_nw - 0.5d0)/areafrac_nw * (0.5d0*dy)
                elseif (areafrac_se < 0.5d0) then  ! CF in cell (i+1,j-1)
                   cf_location(1,axis) = x0(i)   + (0.5d0 - areafrac_se)/(1.0d0 - areafrac_se) * (0.5d0*dx)
                   cf_location(2,axis) = y0(j-1) - (0.5d0 - areafrac_se)/(1.0d0 - areafrac_se) * (0.5d0*dy)
                else   ! CF in cell (i,j)
                   if (areafrac(i,j) >= 0.5d0) then   ! CF in upper left of cell
                      cf_location(1,axis) = x1(i) - (areafrac(i,j) - 0.5d0)/(areafrac(i,j) - areafrac_nw) * (0.5d0*dx)
                      cf_location(2,axis) = y1(j) + (areafrac(i,j) - 0.5d0)/(areafrac(i,j) - areafrac_nw) * (0.5d0*dy)
                   else  ! areafrac(i,j) < 0.5; CF in lower right of cell
                      cf_location(1,axis) = x1(i) + (0.5d0 - areafrac(i,j))/(areafrac_se - areafrac(i,j)) * (0.5d0*dx)
                      cf_location(2,axis) = y1(j) - (0.5d0 - areafrac(i,j))/(areafrac_se - areafrac(i,j)) * (0.5d0*dy)
                   endif
                endif
             endif
          endif
       enddo   ! i
    enddo   ! j

    ! If this proc has a negative value of x, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_location(1,axis), xout=cf_location_ymax, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    if (verbose_calving .and. main_task) then
       print*, ' '
       print*, 'Circular domain: axis, CF location, radius (km)'
       do axis = 1, 8
          radius = sqrt(cf_location(1,axis)**2 + cf_location(2,axis)**2)
          write(6,'(i4,3f10.3)') axis, cf_location(:,axis)/1000.d0, radius/1000.d0
       enddo
    endif

  end subroutine locate_calving_front_circular

!---------------------------------------------------------------------------

  subroutine locate_calving_front_thule(&
       nx,             ny,           &
       dx,             dy,           &
       x0,             y0,           &
       x1,             y1,           &
       parallel,                     &
       itest, jtest,   rtest,        &
       areafrac,                     &
       cf_location)

    use cism_parallel, only: parallel_reduce_maxloc, parallel_reduce_minloc, broadcast
    use glissade_grid_operators, only: glissade_stagger
    use glide_diagnostics, only: point_diag

    ! Find the calving front location along eight profiles on the Thule domain.
    ! These profiles are defined as follows:
    !
    ! Halbrane profiles:
    ! A: (-150,0) to (-150, 740)
    ! B: (150, 0) to ( 150, 740)
    ! C: (-150,0) to (-150,-740)
    ! D: (150, 0) to ( 150,-740)
    !
    ! Caprona profiles:
    ! A: (-390,0) to (-590, 450)
    ! B:  (390,0) to ( 590, 450)
    ! C: (-390,0) to (-590,-450)
    ! D:  (390,0) to ( 590,-450)
    !
    ! The Halbrane profiles are easier. Not sure if I'm going to implement the Caprona profiles.

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy                    ! grid cell size (m)

    real(dp), dimension(nx-1), intent(in) :: x0  ! x coordinate of NE cell corners
    real(dp), dimension(ny-1), intent(in) :: y0  ! y coordinate of NE cell corners
    real(dp), dimension(nx), intent(in) :: x1  ! x coordinate of cell centers
    real(dp), dimension(ny), intent(in) :: y1  ! y coordinate of cell centers

    type(parallel_type), intent(in) :: &
         parallel                  ! info for parallel communication

    real(dp), dimension(nx,ny), intent(in) :: areafrac

    real(dp), dimension(2,8), intent(out) :: &
         cf_location      ! x and y locations of CF along the Halbrane and Caprona profiles
                          ! first index: x and y
                          ! second index: 1 to 4 for Halbrane, 5 to 8 for Caprona

    ! local variables

    integer :: i, j, jj, iglobal, jglobal
    integer :: axis
    integer :: procnum
    real(dp) :: cf_location_xmax, cf_location_ymax, cf_location_xmin, cf_location_ymin, radius
    real(dp) :: this_areafrac, next_areafrac
    real(dp) :: &
         this_areafrac_avg, next_areafrac_avg    ! average of areafrac in two adjacent cells

    real(dp), dimension(nx) :: y_int, areafrac_int
    real(dp) :: x_intercept, y_intercept, slope  ! properties of the profile
    real(dp) :: x_lim, y_lim                     ! outer limits of the profile
    real(dp) :: dist_y, frac_dist

    ! Find the x and y coordinates of the calving front for the Thule domain
    ! along the different profiles specified in CalvingMIP.

    cf_location(:,:) = 0.0d0

    ! Find the CF location along Halbrane profiles A, B, C and D.
    ! All loops are over locally owned cells.
    ! Assume that the x value of each profile coincides with a cell edge
    ! (not a cell center).

    axis = 1  ! index for Halbrane A
    x_intercept = -150.d3

    do i = nhalo+1, nx-nhalo
       if (x0(i) == x_intercept) then  ! E edge of cell lies on the vertical Halbrane profile
          cf_location(1,axis) = x_intercept
          do j = nhalo+1, ny-nhalo
             this_areafrac_avg = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
             next_areafrac_avg = 0.5d0 * (areafrac(i,j+1) + areafrac(i+1,j+1))
             if (this_areafrac_avg > 0.0d0 .and. next_areafrac_avg == 0.0d0) then
                cf_location(2,axis) = y1(j) + (this_areafrac_avg - 0.5d0)*dy
             endif
          enddo
       endif
    enddo

    ! If this proc has a positive value of y, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_location(2,axis), xout=cf_location_ymax, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 2  ! index for Halbrane B (same as A except for positive x_intercept)
    x_intercept = 150.d3

    do i = nhalo+1, nx-nhalo
       if (x0(i) == x_intercept) then  ! E edge of cell lies on the vertical Halbrane profile
          cf_location(1,axis) = x_intercept
          do j = nhalo+1, ny-nhalo
             this_areafrac_avg = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
             next_areafrac_avg = 0.5d0 * (areafrac(i,j+1) + areafrac(i+1,j+1))
             if (this_areafrac_avg > 0.0d0 .and. next_areafrac_avg == 0.0d0) then
                cf_location(2,axis) = y1(j) + (this_areafrac_avg - 0.5d0)*dy
             endif
          enddo
       endif
    enddo

    ! If this proc has a positive value of y, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_location(2,axis), xout=cf_location_ymax, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 3  ! index for Halbrane C (same as A except in the negative y direction)
    x_intercept = -150.d3

    do i = nhalo+1, nx-nhalo
       if (x0(i) == x_intercept) then  ! E edge of cell lies on the vertical Halbrane profile
          cf_location(1,axis) = x_intercept
          do j = ny-nhalo, nhalo+1, -1
             this_areafrac_avg = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
             next_areafrac_avg = 0.5d0 * (areafrac(i,j-1) + areafrac(i+1,j-1))
             if (this_areafrac_avg > 0.0d0 .and. next_areafrac_avg == 0.0d0) then
                cf_location(2,axis) = y1(j) + (0.5d0 - this_areafrac_avg)*dy
             endif
          enddo
       endif
    enddo

    ! If this proc has a negative value of y, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_location(2,axis), xout=cf_location_ymin, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 4  ! index for Halbrane D (same as C except for positive x_intercept)
    x_intercept = 150.d3

    do i = nhalo+1, nx-nhalo
       if (x0(i) == x_intercept) then  ! E edge of cell lies on the vertical Halbrane profile
          cf_location(1,axis) = x_intercept
          do j = ny-nhalo, nhalo+1, -1
             this_areafrac_avg = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
             next_areafrac_avg = 0.5d0 * (areafrac(i,j-1) + areafrac(i+1,j-1))
             if (this_areafrac_avg > 0.0d0 .and. next_areafrac_avg == 0.0d0) then
                cf_location(2,axis) = y1(j) + (0.5d0 - this_areafrac_avg)*dy
             endif
          enddo
       endif
    enddo

    ! If this proc has a negative value of y, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_location(2,axis), xout=cf_location_ymin, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    ! Find the CF location along Caprona profiles A, B, C and D.
    ! The Caprona profiles cut across cells without passing through centers or corners.
    ! As a result, the logic below is more complicated than for the Halbrane profiles,
    !  and more approximate. Results might be better with offline interpolation of areafrac.

    axis = 5  ! index for Caprona A
    x_intercept = -390.d3
    x_lim = -590.d3
    y_lim = 450.d3
    slope = y_lim/(x_lim - x_intercept)  ! rise over run = 450/(-200) = -2.25
    y_intercept = -x_intercept * slope

    ! Adjust x_lim to allow the CF to be a little out of bounds
    x_lim = -650.d3

    y_int = 0.0d0
    areafrac_int = 0.0d0

    ! Estimate areafrac at each point where the Caprona profile intersects the x1 grid
    do i = nx-nhalo, nhalo+1, -1
       if (x1(i) < x_intercept .and. x1(i) >= x_lim) then  ! x1 in range
          y_int(i) = slope*x1(i) + y_intercept  ! profile intersects x1 grid at (x1(i),y)
          do j = nhalo+1, ny-nhalo
             if (y_int(i) >= y1(j) .and. y_int(i) < y1(j+1)) then
                ! Interpolate to estimate a_eff at (x1(i),y_int)
                areafrac_int(i) = areafrac(i,j) + (y_int(i) - y1(j))/dy * (areafrac(i,j+1) - areafrac(i,j))
                exit
             endif
          enddo
       endif
    enddo

    ! Find a point along the profile where the interpolated areafrac = 0.5
    do i = nx-nhalo, nhalo+1, -1
       if (areafrac_int(i) > 0.5d0 .and. areafrac_int(i-1) < 0.5d0) then
          dist_y = y_int(i-1) - y_int(i) ! y distance between neighboring intersection points
          frac_dist = (areafrac_int(i) - 0.5d0) / (areafrac_int(i) - areafrac_int(i-1))
          cf_location(1,axis) = x1(i) - frac_dist*dx
          cf_location(2,axis) = y_int(i) + frac_dist*dist_y
          exit
       endif
    enddo

    ! If this proc has a positive value of y, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_location(2,axis), xout=cf_location_ymin, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 6  ! index for Caprona B
    x_intercept = 390.d3
    x_lim = 590.d3
    y_lim = 450.d3
    slope = y_lim/(x_lim - x_intercept)  ! rise over run = 450/200 = 2.25
    y_intercept = -x_intercept * slope

    ! Adjust x_lim to allow the CF to be a little out of bounds
    x_lim = 650.d3

    y_int = 0.0d0
    areafrac_int = 0.0d0

    ! Estimate areafrac at each point where the Caprona profile intersects the x1 grid
    do i = nhalo+1, nx-nhalo
       if (x1(i) >= x_intercept .and. x1(i) < x_lim) then  ! x1 in range
          y_int(i) = slope*x1(i) + y_intercept  ! profile intersects x1 grid at (x1(i),y)
          do j = nhalo+1, ny-nhalo
             if (y_int(i) >= y1(j) .and. y_int(i) < y1(j+1)) then
                ! Interpolate to estimate a_eff at (x1(i),y_int)
                areafrac_int(i) = areafrac(i,j) + (y_int(i) - y1(j))/dy * (areafrac(i,j+1) - areafrac(i,j))
                exit
             endif
          enddo
       endif
    enddo

    if (verbose_calving .and. this_rank ==rtest) then
!       print*, 'Caprona B intersection points: i, x, y, areafrac'
!       do i = nhalo+1, nx-nhalo
!          if (y_int(i) /= 0.0d0) then
!             write(6,'(i4,3f10.3)'), i, x1(i)/1000.d0, y_int(i)/1000.d0, areafrac_int(i)
!          endif
!       enddo
    endif

    ! Find a point along the profile where the interpolated areafrac = 0.5
    do i = nhalo+1, nx-nhalo
       if (areafrac_int(i) > 0.5d0 .and. areafrac_int(i+1) < 0.5d0) then
          dist_y = y_int(i+1) - y_int(i) ! y distance between neighboring intersection points
          frac_dist = (areafrac_int(i) - 0.5d0) / (areafrac_int(i) - areafrac_int(i+1))
          cf_location(1,axis) = x1(i) + frac_dist*dx
          cf_location(2,axis) = y_int(i) + frac_dist*dist_y
          if (verbose_calving .and. this_rank == rtest) then
!             print*, '1st IP: x, y, a_eff =', x1(i)/1000.d0, y_int(i)/1000.d0, areafrac_int(i)
!             print*, '2nd IP: x, y, a_eff =', x1(i+1)/1000.d0, y_int(i+1)/1000.d0, areafrac_int(i+1)
!             print*, 'dist_y, frac_dist =', dist_y/1000.d0, frac_dist
!             print*, 'CF location =', cf_location(1,axis)/1000.d0, cf_location(2,axis)/1000.d0
!             print*, 'residual y - (mx + b):', cf_location(2,axis) - slope*cf_location(1,axis) - y_intercept
          endif
          exit
       endif
    enddo

    ! If this proc has a positive value of y, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_location(2,axis), xout=cf_location_ymin, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 7  ! index for Caprona C
    x_intercept = -390.d3
    x_lim = -590.d3
    y_lim = -450.d3
    slope = y_lim/(x_lim - x_intercept)  ! rise over run = -450/(-200) = 9/4
    y_intercept = -x_intercept * slope

    ! Adjust x_lim to allow the CF to be a little out of bounds
    x_lim = -650.d3

    y_int = 0.0d0
    areafrac_int = 0.0d0

    ! Estimate areafrac at each point where the Caprona profile intersects the x1 grid
    do i = nx-nhalo, nhalo+1, -1
       if (x1(i) < x_intercept .and. x1(i) >= x_lim) then  ! x1 in range
          y_int(i) = slope*x1(i) + y_intercept  ! profile intersects x1 grid at (x1(i),y)
          do j = ny-nhalo, nhalo+1, -1
             if (y_int(i) <= y1(j) .and. y_int(i) > y1(j-1)) then
                ! Interpolate to estimate a_eff at (x1(i),y_int)
                areafrac_int(i) = areafrac(i,j) + (y1(j) - y_int(i))/dy * (areafrac(i,j-1) - areafrac(i,j))
                exit
             endif
          enddo
       endif
    enddo

    ! Find a point along the profile where the interpolated areafrac = 0.5
    do i = nx-nhalo, nhalo+1, -1
       if (areafrac_int(i) > 0.5d0 .and. areafrac_int(i-1) < 0.5d0) then
          dist_y = y_int(i-1) - y_int(i) ! y distance between neighboring intersection points
          frac_dist = (areafrac_int(i) - 0.5d0) / (areafrac_int(i) - areafrac_int(i-1))
          cf_location(1,axis) = x1(i) - frac_dist*dx
          cf_location(2,axis) = y_int(i) + frac_dist*dist_y
          exit
       endif
    enddo

    ! If this proc has a negative value of y, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_location(2,axis), xout=cf_location_ymin, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    axis = 8  ! index for Caprona D
    x_intercept = 390.d3
    x_lim = 590.d3
    y_lim = -450.d3
    slope = y_lim/(x_lim - x_intercept)  ! rise over run = -450/200 = -9/4
    y_intercept = -x_intercept * slope

    ! Adjust x_lim to allow the CF to be a little out of bounds
    x_lim = 650.d3

    y_int = 0.0d0
    areafrac_int = 0.0d0

    ! Estimate areafrac at each point where the Caprona profile intersects the x1 grid
    do i = nhalo+1, nx-nhalo
       if (x1(i) >= x_intercept .and. x1(i) < x_lim) then  ! x1 in range
          y_int(i) = slope*x1(i) + y_intercept  ! profile intersects x1 grid at (x1(i),y)
          do j = ny-nhalo, nhalo+1, -1
             if (y_int(i) <= y1(j) .and. y_int(i) > y1(j-1)) then
                ! Interpolate to estimate a_eff at (x1(i),y_int)
                areafrac_int(i) = areafrac(i,j) + (y1(j) - y_int(i))/dy * (areafrac(i,j-1) - areafrac(i,j))
                exit
             endif
          enddo
       endif
    enddo

    ! Find a point along the profile where the interpolated areafrac = 0.5
    do i = nhalo+1, nx-nhalo
       if (areafrac_int(i) > 0.5d0 .and. areafrac_int(i+1) < 0.5d0) then
          dist_y = y_int(i+1) - y_int(i) ! y distance between neighboring intersection points
          frac_dist = (areafrac_int(i) - 0.5d0) / (areafrac_int(i) - areafrac_int(i+1))
          cf_location(1,axis) = x1(i) + frac_dist*dx
          cf_location(2,axis) = y_int(i) + frac_dist*dist_y
       endif
    enddo

    ! If this proc has a negative value of y, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_location(2,axis), xout=cf_location_ymin, xprocout=procnum)
    call broadcast(cf_location(:,axis), proc=procnum)

    if (verbose_calving .and. main_task) then
       print*, ' '
       print*, 'Thule domain: axis, CF location, radius (km)'
       do axis = 1, 8
          radius = sqrt(cf_location(1,axis)**2 + cf_location(2,axis)**2)
          write(6,'(i4,3f10.3)') axis, cf_location(:,axis)/1000.d0, radius/1000.d0
       enddo
    endif

  end subroutine locate_calving_front_thule

!---------------------------------------------------------------------------

end module glissade_calving

!---------------------------------------------------------------------------
