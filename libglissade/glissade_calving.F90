!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_calving.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2017
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
  use parallel

  !WHL - debug
  use glimmer_paramets, only: thk0

  implicit none

  ! colors for fill subroutine
  integer, parameter :: initial_color = 0   ! initial color, represented by integer
  integer, parameter :: fill_color = 1      ! fill color, represented by integer
  integer, parameter :: boundary_color = -1 ! boundary color, represented by integer

  !WHL - debug
!!  logical, parameter :: verbose_calving = .false.
  logical, parameter :: verbose_calving = .true.

contains

!-------------------------------------------------------------------------------  

  subroutine glissade_calving_mask_init(which_calving,                       &
                                        dx,                dy,               &
                                        calving_front_x,   calving_front_y,  &
                                        calving_mask)

    ! Compute an integer calving mask if needed for the CALVING_GRID_MASK option

    ! Input/output arguments

    integer,  intent(in) :: which_calving       !> option for calving law
    real(dp), intent(in) :: dx, dy              !> cell dimensions in x and y directions (m)
    real(dp), intent(in) :: calving_front_x     !> for CALVING_GRID_MASK option, calve ice wherever abs(x) > calving_front_x (m)
    real(dp), intent(in) :: calving_front_y     !> for CALVING_GRID_MASK option, calve ice wherever abs(y) > calving_front_y (m)
    integer, dimension(:,:), intent(inout) :: calving_mask   !> output mask: calve floating ice wherever the mask = 1

    ! Local variables

    real(dp) :: xcell, ycell     ! global cell center coordinates (m)
    integer :: nx, ny            ! horizontal grid dimensions
    integer :: i, j              ! local cell indices
    integer :: iglobal, jglobal  ! global cell indices

    nx = size(calving_mask,1)
    ny = size(calving_mask,2)

    if (which_calving == CALVING_GRID_MASK) then  ! calve where calving_mask = 1

       ! Two possibilities to consider:
       ! (1) The calving mask was read from an input file.
       !     If so, then there should be at least one grid location where calving_mask = 1.
       !     There is nothing more to do here.
       ! (2) The calving mask was not read from an input file, but is determined here based on the values
       !     of calving_front_x and calving_front_y.

       if (maxval(calving_mask) > 0) then

          ! calving_mask was read from the input file; do not need to compute a mask here

       else   ! create calving_mask here

          ! initialize
          calving_mask(:,:) = 0   ! no calving by default

          if (calving_front_x > 0.0d0) then

             ! set calving_mask = 1 where abs(x) > calving_front_x

             do j = 1, ny
                do i = 1, nx

                   ! find global i and j indices
                   call parallel_globalindex(i, j, iglobal, jglobal)

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
                   call parallel_globalindex(i, j, iglobal, jglobal)

                   ! find cell center y coordinate
                   ycell = (dble(jglobal) - 0.5d0) * dy

                   ! set calving mask = 1 based on cell coordinates relative to the calving front
                   if (abs(ycell) > calving_front_y) then
                      calving_mask(i,j) = 1
                   endif

                enddo   ! i
             enddo   ! j

          endif   ! calving_front_y > 0

       endif   ! maxval(calving_mask) > 0

    else  ! whichcalving /= CALVING_GRID_MASK

       ! do nothing; other calving options do not require initialization

    endif

  end subroutine glissade_calving_mask_init

!-------------------------------------------------------------------------------

  subroutine glissade_calve_ice(which_calving,     calving_domain,   &
                                itest,    jtest,   rtest,            &
                                thck,              relx,             &
                                topg,              eus,              &
                                thklim,                  &
                                marine_limit,            &
                                calving_fraction,        &
                                calving_timescale,       &
                                dt,                      &
                                efvs_vertavg,            &
                                strain_rate_determinant, &
                                eigencalving_constant,   &
                                calving_minthck_constant,&
                                calving_mask,            &
                                remove_floating_islands, &
                                floating_path_minthck,   &
                                damage,                  &
                                damage_threshold,        &
                                damage_column,           &
                                sigma,                   &
                                calving_thck)

    ! Calve ice according to one of several methods

    use glissade_masks

    implicit none

    !---------------------------------------------------------------------
    ! Subroutine arguments
    ! Currently, thck, relx, topg, eus, marine_limit, calving_minthck and calving_thck are scaled by thk0
    ! TODO - Remove thickness scaling from this subroutine.  Would need to multiply several input arguments by thk0. 
    !---------------------------------------------------------------------

    integer,  intent(in)  :: which_calving        !> option for calving law
    integer,  intent(in)  :: calving_domain       !> option for where calving can occur
                                                  !> = 0 if calving occurs at the ocean edge only
                                                  !> = 1 if calving occurs everywhere the calving criterion is met
                                                  !> = 2 if calving occurs where criterion is met and there is a connected path
                                                  !>     to the ocean through other cells where the criterion is met
    integer, intent(in)   :: itest, jtest, rtest  !> coordinates of diagnostic point

    real(dp), dimension(:,:), intent(inout) :: thck              !> ice thickness
    real(dp), dimension(:,:), intent(in)    :: relx              !> relaxed bedrock topography
    real(dp), dimension(:,:), intent(in)    :: topg              !> present bedrock topography
    real(dp), intent(in)                    :: eus               !> eustatic sea level
    real(dp), intent(in)                    :: thklim            !> minimum thickness for dynamically active ice; used by remove_floating_islands
    real(dp), intent(in)                    :: marine_limit      !> lower limit on topography elevation at marine edge before ice calves
    real(dp), intent(in)                    :: calving_fraction  !> fraction of ice lost at marine edge when calving; 
                                                                 !> used with which_ho_calving = CALVING_FLOAT_FRACTION
    real(dp), intent(in)                    :: calving_timescale !> time scale for calving; calving_thck = thck * max(dt/calving_timescale, 1)
                                                                 !> if calving_timescale = 0, then calving_thck = thck
    real(dp), intent(in)                    :: dt                !> model timestep (used with calving_timescale)
    real(dp), dimension(:,:), intent(in)    :: efvs_vertavg      !> vertical mean effective viscosity (Pa s); used for eigencalving
    real(dp), dimension(:,:), intent(in)    :: strain_rate_determinant  !> determinant of horizontal strain rate tensor (s^{-2})
                                                                        !> used for eigencalving
    real(dp), intent(in)                    :: eigencalving_constant    !> dimensionless constant used for eigencalving
    real(dp), intent(in)                    :: calving_minthck_constant !> min thickness of floating ice before it calves;
                                                                        !> used with which_ho_calving = CALVING_THCK_THRESHOLD
    integer,  dimension(:,:), intent(in)    :: calving_mask            !> integer mask: calve ice where mask = 1
    logical,  intent(in)                    :: remove_floating_islands !> if true, then remove floating ice islands
    real(dp), intent(in)                    :: floating_path_minthck   !> min thickness along path from floating ice back to grounded ice
!    real(dp), dimension(:,:,:), intent(in)  :: damage            !> 3D scalar damage parameter
    real(dp), dimension(:,:,:), intent(inout)  :: damage         !> 3D scalar damage parameter  !WHL - 'inout' if damage is updated below
    real(dp), dimension(:,:), intent(out)   :: damage_column     !> 2D vertically integrated scalar damage parameter
    real(dp), intent(in)                    :: damage_threshold  !> threshold value where ice is sufficiently damaged to calve
    real(dp), dimension(:), intent(in)      :: sigma             !> vertical sigma coordinate
    real(dp), dimension(:,:), intent(out)   :: calving_thck      !> thickness lost due to calving in each grid cell

    integer :: nx, ny      ! horizontal grid dimensions
    integer :: nz          ! number of vertical levels
                           ! Note: number of ice layers = nz-1
    integer :: i, j, k
    integer :: count, maxcount_fill  ! loop counters

    integer,  dimension(:,:), allocatable   ::  &
         color                ! integer 'color' for filling the calving domain (with CALVING_DOMAIN_OCEAN_CONNECT)

    ! basic masks
    integer, dimension(:,:), allocatable   ::  &
         ice_mask,               & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,          & ! = 1 where ice is present and floating, else = 0
         ocean_mask                ! = 1 where topg is below sea level and thk <= thklim, else = 0

    ! masks for CALVING_THCK_THRESHOLD
    integer, dimension(:,:), allocatable   ::  &
         thick_or_grounded_mask, & ! = 1 for grounded ice and thick floating ice, else = 0
         protected_floating_mask   ! = 1 for thin floating ice in protected ring, else = 0

    ! other calving-specific fields
    real(dp), dimension(:,:), allocatable  ::  &
         calving_minthck      ! min thickness (m) of floating ice before it calves

    ! Note: Calving occurs in a cell if and only if (1) the calving law permits calving, 
    !       and (2) the cell is in the calving domain, as specified by the calving_domain option.
    !       The calving domain by default is limited to the ocean edge (CALVING_DOMAIN_OCEAN_EDGE), 
    !       but can be extended to include all ice-covered cells (CALVING_DOMAIN_EVERYWHERE), or
    !       cells connected to the ocean through other cells that meet the calving criterion
    !       (CALVING_DOMAIN_OCEAN_CONNECT).
    ! TODO: Change the default to calving_domain_everywhere?

    !TODO - Make these integer masks like the ones above?
    logical, dimension(:,:), allocatable   ::  &
         calving_law_mask,    & ! = T where the calving law permits calving, else = F
         calving_domain_mask    ! = T in the domain where calving is allowed to occur (e.g., at ocean edge), else = F

    real(dp) :: &
         float_fraction_calve  ! = calving_fraction for which_calving = CALVING_FLOAT_FRACTION
                               ! = 1.0 for which_calving = CALVING_FLOAT_ZERO
   
    real(dp) :: &
         tau_spreading         ! stress associated with horizontal spreading (Pa)

    !WHL - debug
    integer :: sum_fill_local, sum_fill_global  ! number of filled cells

    integer :: iplot1, iplot2

    !default
!    iplot1 = nx-20
!    iplot2 = nx-1

    ! initialize
    calving_thck(:,:) = 0.d0

    nx = size(thck,1)
    ny = size(thck,2)
    nz = size(sigma)


    if (which_calving == CALVING_NONE) then

       ! remove floating islands (if desired) and return

       if (remove_floating_islands) then

          if (verbose_calving .and. main_task) then
             print*, 'Remove floating islands'
          endif

          ! Set the thickness thresholds for cells to be counted as active.
          ! Floating regions will be identified as ice islands unless they are connected to grounded ice
          !  along a path consisting only of active cells (i.e., cells with thk > thklim).
          ! Note: A limit of 0.0 does not work because it erroneously counts very thin floating cells as active.
          !       Then the algorithm can fail to identify floating regions that should be removed.

          call glissade_remove_floating_islands(&
               thck,          relx,                  &
               topg,          eus,                   &
               thklim,        floating_path_minthck, &
               calving_thck)

       endif

       return

    endif

    allocate (ice_mask(nx,ny))
    allocate (floating_mask(nx,ny))
    allocate (ocean_mask(nx,ny))

    allocate (calving_law_mask(nx,ny))
    allocate (calving_domain_mask(nx,ny))

    !WHL - debug
    if (verbose_calving .and. main_task) then
       print*, ' '
       print*, 'In glissade_calve_ice'
       print*, 'which_calving =', which_calving
       print*, 'calving_domain =', calving_domain
!       print*, 'i, relx, topg, thck, usfc:'
!       do i = iplot1, iplot2
!          print*, i, relx(i,j)*thk0, topg(i,j)*thk0, thck(i,j)*thk0, (topg(i,j) + thck(i,j))*thk0
!       enddo
    endif

    ! Set the thickness fraction to be removed in each calving cell
    ! Note: The CALVING_FLOAT_FRACTION option has been superseded by the calving_timescale variable,
    !       but is included here for consistency with Glide.

    if (which_calving == CALVING_FLOAT_FRACTION) then

       !WHL - Changed definition of calving fraction; now it is the fraction lost
       !      rather than the fraction remaining
       float_fraction_calve = calving_fraction
       
    else  ! other calving options

       if (calving_timescale == 0.0d0) then  ! calve the entire column for eligible columns (this is the default)
          float_fraction_calve = 1.0d0
       else  ! calve a fraction of the column based on the calving time scale
          float_fraction_calve = min(dt/calving_timescale, 1.0d0)
       endif
       
    endif
       
    ! Get masks
    ! Use thickness limit of 0.0 instead of thklim so as to remove ice from any cell
    !  that meets the calving criteria, not just dynamically active ice

    call glissade_get_masks(nx,            ny,            &
                            thck,          topg,          &
                            eus,           0.0d0,         &   ! thklim = 0.0
                            ice_mask,      floating_mask, &
                            ocean_mask)

    if (verbose_calving .and. this_rank==rtest) then

       print*, ' '
       print*, 'thck field, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+2, jtest-2, -1
          write(6,'(i6)',advance='no') j
          do i = itest-2, itest+2
             write(6,'(f10.3)',advance='no') thck(i,j)*thk0
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'ice_mask, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+2, jtest-2, -1
          write(6,'(i6)',advance='no') j
          do i = itest-2, itest+2
             write(6,'(L10)',advance='no') ice_mask(i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'floating_mask, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+2, jtest-2, -1
          write(6,'(i6)',advance='no') j
          do i = itest-2, itest+2
             write(6,'(L10)',advance='no') floating_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '

    endif

    ! For thickness-based calving (including eigencalving), compute the 2D field calving_minthck.
    ! Ice thinner than calving_minthck satisfies the calving criterion.

    allocate(calving_minthck(nx,ny))
    calving_minthck(:,:) = 0.0d0

    if (which_calving == EIGENCALVING) then

       ! Compute the minimum thickness using an eigencalving criterion.
       ! The assumption is that ice calves when
       !
       ! k * 2*efvs*sqrt(determinant) > 0.5*rhoi*grav*thck
       !
       ! where efvs = vertical mean effective viscosity
       !       determinant = determinant of the horizontal strain rate tensor
       !       k ~ 10 is a dimensionless empirical constant accounting for damage
       !
       ! The LHS has units of stress; it is positive only when both eigenvalues are positive
       !  (i.e., the ice is spreading in both principal directions).
       ! The RHS is the mean value of the hydrostatic stress.
       ! This calving law is based on the eigencalving scheme of Levermann et al. (2012).
       ! The difference is that instead of computing a calving rate, we compute a thickness threshold
       !  based on the assumption that thinner ice with lower hydrostatic pressure calves more easily.
       ! The determinant of the strain rate tensor is computed after the velocity calculation and
       !  has units of s^{-2}.
       
       do j = 1, ny
          do i = 1, nx
             
             if (strain_rate_determinant(i,j) > 0.0d0) then
                tau_spreading = 2.d0 * efvs_vertavg(i,j) * sqrt(strain_rate_determinant(i,j))   ! Pa
                calving_minthck(i,j) = eigencalving_constant * tau_spreading / (0.5d0 * rhoi * grav)
             else
                calving_minthck(i,j) = 0.0d0
             endif
             
             ! calving_minthck has units of m; scale to model units
             calving_minthck(i,j) = calving_minthck(i,j) / thk0

          enddo  ! j
       enddo  ! i

       !WHL - Debug
       if (verbose_calving .and. this_rank == rtest) then

          print*, ' '
          print*, 'Eigencalving fields:'

          print*, ' '
          print*, 'sqrt(determinant) (yr-1), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+2, jtest-2, -1
             write(6,'(i6)',advance='no') j
             do i = itest-2, itest+2
                write(6,'(e10.3)',advance='no') sqrt(strain_rate_determinant(i,j)) * scyr
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'efvs_vertavg (Pa yr), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+2, jtest-2, -1
             write(6,'(i6)',advance='no') j
             do i = itest-2, itest+2
                write(6,'(e10.3)',advance='no') efvs_vertavg(i,j) * scyr
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'tau_spreading (Pa), itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+2, jtest-2, -1
             write(6,'(i6)',advance='no') j
             do i = itest-2, itest+2
                write(6,'(e10.3)',advance='no') 2.d0 * efvs_vertavg(i,j) * sqrt(strain_rate_determinant(i,j))
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'mean hydrostatic pressure, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+2, jtest-2, -1
             write(6,'(i6)',advance='no') j
             do i = itest-2, itest+2
                write(6,'(e10.3)',advance='no') 0.5d0 * rhoi * grav * thck(i,j) *thk0
             enddo
             write(6,*) ' '
          enddo
          
          print*, ' '
          print*, 'thck, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+2, jtest-2, -1
             write(6,'(i6)',advance='no') j
             do i = itest-2, itest+2
                write(6,'(f10.3)',advance='no') thck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo

          print*, ' '
          print*, 'calving_minthck, itest, jtest, rank =', itest, jtest, rtest
          do j = jtest+2, jtest-2, -1
             write(6,'(i6)',advance='no') j
             do i = itest-2, itest+2
                write(6,'(f10.3)',advance='no') calving_minthck(i,j)*thk0
             enddo
             write(6,*) ' '
          enddo
          
       endif
       
    elseif (which_calving == CALVING_THCK_THRESHOLD) then
       
       ! set to a constant everywhere
       calving_minthck(:,:) = calving_minthck_constant
       
    endif  ! thickness-based calving

    ! Do the calving based on the value of which_calving

    ! Note: The thickness-threshold option is different from the others.
    !       For the other options, we look at each cell and determine whether it meets the calving-law criteria
    !        (e.g., ice is floating, or the topography lies below a given level).
    !        If a cell meets the criteria and lies in the calving domain (e.g., at the margin), it is calved.
    !       For the thickness-threshold option, ice thinner than calving_minthck is calved, but only if it 
    !        lies beyond a protected ring of thin floating ice at the margin.
    !       The reason for this more complicated approach is that we do not want to remove all floating ice
    !         with thck < calving_minthck, because then we would remove thin ice that has just been advected
    !         from active cells at the margin, and thus the calving front would not be able to advance.
    !       By protecting a ring of inactive ice (thck < calving_minthck) at the margin, we allow ice in
    !        these cells to thicken, thus advancing the calving front.
    !       The calving front can retreat when floating ice at the margin thins to a value less than calving_minthck,
    !        removing protection from previously protected cells.

    if (which_calving == CALVING_THCK_THRESHOLD .or. which_calving == EIGENCALVING) then
       ! calve floating ice thinner than calving_minthck (if not part of the protected ring)

       ! initialize masks

       allocate (thick_or_grounded_mask(nx,ny))
       allocate (protected_floating_mask(nx,ny))

       thick_or_grounded_mask(:,:) = 0
       protected_floating_mask(:,:) = 0

       ! Identify cells that are (1) grounded or (2) floating with thck > calving_minthck.
       ! These cells are protected from calving.
       do j = 1, ny
          do i = 1, nx
             if (ice_mask(i,j) == 1) then
!!                if (floating_mask(i,j) == 0 .or. thck(i,j) > calving_minthck) then
                if (floating_mask(i,j) == 0 .or. thck(i,j) > calving_minthck(i,j)) then
                   thick_or_grounded_mask(i,j) = 1
                endif
             endif
          enddo  ! i
       enddo  ! j

       ! Identify cells in the protected ring of thin floating ice.
       ! Thin floating cells are protected if adjacent to cells that were protected above.
       do j = 2, ny-1
          do i = 2, nx-1
!!             if (floating_mask(i,j) == 1 .and. thck(i,j) <= calving_minthck) then
             if (floating_mask(i,j) == 1 .and. thck(i,j) <= calving_minthck(i,j)) then
                if (thick_or_grounded_mask(i-1,j) == 1 .or. thick_or_grounded_mask(i+1,j) == 1 .or.  &
                    thick_or_grounded_mask(i,j-1) == 1 .or. thick_or_grounded_mask(i,j+1) == 1) then
                   protected_floating_mask(i,j) = 1
                endif
             endif
          enddo  ! i
       enddo  ! j

       ! halo update, since the loop above misses some halo cells
       call parallel_halo(protected_floating_mask)

       ! Calve thin floating ice, except in the protected ring.
       ! Note: The calving_law_mask does not need to be set because it is not used subsequently,
       !       but could be set here for diagnostic purposes.

       do j = 1, ny
          do i = 1, nx
             if (floating_mask(i,j) == 1  .and. thck(i,j) <= calving_minthck(i,j) .and. protected_floating_mask(i,j) == 0) then
                if (verbose_calving .and. thck(i,j) > 0.0d0) &
                     print*, 'Calve thin floating ice: task, i, j, thck =', this_rank, i, j, thck(i,j)*thk0
                calving_thck(i,j) = float_fraction_calve * thck(i,j)
                thck(i,j) = thck(i,j) - calving_thck(i,j)
                !WHL - Also handle tracers?  E.g., set damage(:,i,j) = 0.d0?
             endif
          enddo
       enddo

    else   ! other calving options

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

          where (relx <= marine_limit + eus)
             calving_law_mask = .true.
          elsewhere
             calving_law_mask = .false.
          endwhere

       case(CALVING_TOPG_THRESHOLD)   ! set thickness to zero if present bedrock is below a given level

          where (topg < marine_limit + eus)
             calving_law_mask = .true.
          elsewhere
             calving_law_mask = .false.
          endwhere

       case(CALVING_GRID_MASK)     ! calve ice where the input calving mask = 1

          where (calving_mask == 1)
             calving_law_mask = .true.
          elsewhere
             calving_law_mask = .false.
          endwhere

       case(CALVING_HUYBRECHTS)    ! Huybrechts grounding line scheme for Greenland initialization

       !WHL - Previously, this code assumed that eus and relx have units of meters.
       !      Changed to be consistent with dimensionless thickness units.
!       if (eus > -80.d0) then
!          where (relx <= 2.d0*eus)
!             calving_thck = thck
!             thck = 0.0d0
!          end where
!       elseif (eus <= -80.d0) then
!          where (relx <= (2.d0*eus - 0.25d0*(eus + 80.d0)**2.d0))
!             calving_thck = thck
!             thck = 0.0d0
!          end where
!       end if
          if (eus*thk0 > -80.d0) then
             where (relx*thk0 <= 2.d0*eus*thk0)
                calving_law_mask = .true.
             elsewhere
                calving_law_mask = .false.
             end where
          elseif (eus*thk0 <= -80.d0) then
             where (relx*thk0 <= (2.d0*eus*thk0 - 0.25d0*(eus*thk0 + 80.d0)**2.d0))
                calving_law_mask = .true.
             elsewhere
                calving_law_mask = .false.
             end where
          end if

       case(CALVING_DAMAGE)   ! remove ice that is sufficiently damaged
                              !WHL - This is a rough initial implementation

          !WHL - debug - test damage field for MISMIP
!          print*, 'nx, ny =', nx, ny
!          print*, 'Prescribe damage:'
!          damage(:,:,:) = 0.d0
!          do j = nhalo+1, ny-nhalo
!             do i = nhalo+1, nx-nhalo
!                if (j == 3) then
!                   damage(:,i,j) = damage_threshold + 0.1d0
!                endif
!             enddo
!          enddo

          !TODO - Define ice_mask = 1 for dynamically active ice only (thck > thklim)?

          ! Diagnose the vertically integrated damage in each column,
          ! assuming the 3D damage field has been prognosed external to this subroutine.
          !WHL - For now, simply compute the thickness-weighted mean damage
          damage_column(:,:) = 0.0d0
          do j = 1, ny
             do i = 1, nx
                if (ice_mask(i,j) == 1) then
                   do k = 1, nz-1
                      damage_column(i,j) = damage_column(i,j) + damage(k,i,j) * (sigma(k+1) - sigma(k))
                   enddo
                endif
             enddo
          enddo
          
          ! set calving-law mask based on the vertically integrated damage
          where (damage_column > damage_threshold)  !WHL - could use '>=' instead of '>' if preferred
             calving_law_mask = .true.
          elsewhere
             calving_law_mask = .false.
          endwhere

          !WHL - debug - print values of calving_law_mask
!          if (main_task) then
!             print*, 'i, j, damage, calving_law_mask:'
!             j = 3
!             do i = 1, nx
!                print*, i, j, damage_column(i,j), calving_law_mask(i,j)
!             enddo
!          endif

       end select

       ! halo update (may not be necessary if thck, damage, etc. are correct in halos, but including to be safe)
       call parallel_halo(calving_law_mask)

       ! set the calving domain mask

       if (calving_domain == CALVING_DOMAIN_OCEAN_EDGE) then  ! calving domain includes floating cells at margin only
                                                              !WHL - Could modify to include grounded marine cells at margin
          do j = 2, ny-1
             do i = 2, nx-1

                if (verbose_calving .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                   print*, 'task, i, j, ice_mask, floating_mask:',  &
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
          call parallel_halo(calving_domain_mask)

       elseif (calving_domain == CALVING_DOMAIN_EVERYWHERE) then  ! calving domain includes all cells

          calving_domain_mask(:,:) = .true.

       elseif (calving_domain == CALVING_DOMAIN_OCEAN_CONNECT) then

          calving_domain_mask(:,:) = .false.

          ! initialize
          ! Assign the initial color to cells that meet the calving-law criteria and thus could calve,
          !  but only if they are connected to the ocean through other cells that meet the criteria.
          ! Assign the boundary color to cells that do not meet the calving-law criteria.

          allocate (color(nx,ny))
          do j = 1, ny
             do i = 1, nx
                if (calving_law_mask(i,j)) then
                   color(i,j) = initial_color
                else
                   color(i,j) = boundary_color
                endif
             enddo
          enddo

          ! Loop through cells, identifying cells that lie on the ocean margin.
          ! Fill each such cell with calving_law_mask = T, and recursively fill neighbor cells with calving_law_mask = T.
          ! We may have to do this several times to incorporate connections between neighboring processors.
          ! Alternatively, we could use a smaller value of maxcount_fill and allow calving to occur over multiple time steps,
          !  but this could make the calving evolution more dependent on ntasks than desired.

!          maxcount_fill = 1  ! setting maxcount_fill = 1 will ignore connections between neighboring processors.
          maxcount_fill = max(ewtasks,nstasks)

          do count = 1, maxcount_fill

             if (count == 1) then   ! identify margin cells that can seed the fill

                do j = 1, ny
                   do i = 1, nx
                      if ( floating_mask(i,j) == 1 .and. calving_law_mask(i,j) .and.  &
                           (ocean_mask(i-1,j)==1 .or. ocean_mask(i+1,j)==1 .or. ocean_mask(i,j-1)==1 .or. ocean_mask(i,j+1)==1) ) then
                         if (color(i,j) /= boundary_color .and. color(i,j) /= fill_color) then
                            ! assign the fill color to this cell, and recursively fill neighbor cells
                            call glissade_fill(nx, ny, i, j, color)
                         endif
                      endif
                   enddo
                enddo

             else  ! count > 1; check for halo cells that were just filled on neighbor processors

                call parallel_halo(color)

                ! west halo layer
                i = nhalo
                do j = 1, ny
                   if (color(i,j) == fill_color) call glissade_fill(nx, ny, i+1, j, color)
                enddo

                ! east halo layer
                i = nx - nhalo + 1
                do j = 1, ny
                   if (color(i,j) == fill_color) call glissade_fill(nx, ny, i-1, j, color)
                enddo

                ! south halo layer
                j = nhalo
                do i = nhalo+1, nx-nhalo  ! already checked halo corners above
                   if (color(i,j) == fill_color) call glissade_fill(nx, ny, i, j+1, color)
                enddo

                ! north halo layer
                j = ny-nhalo+1
                do i = nhalo+1, nx-nhalo  ! already checked halo corners above
                   if (color(i,j) == fill_color) call glissade_fill(nx, ny, i, j-1, color)
                enddo

             endif  ! count = 1

             sum_fill_local = 0
             do j = nhalo+1, ny-nhalo
                do i = nhalo+1, nx-nhalo
                   if (color(i,j) == fill_color) sum_fill_local = sum_fill_local + 1
                enddo
             enddo

             !WHL - If running a large problem, may want to reduce the frequency of this global sum
             sum_fill_global = parallel_reduce_sum(sum_fill_local)

             if (verbose_calving) then
                print*, 'this_rank, sum_fill_local, sum_fill_global:', this_rank, sum_fill_local, sum_fill_global 
             endif

          enddo  ! count

          ! At this point, all cells with calving_law_mask = T should have the fill color if they
          !  are connected to the margin through other cells with calving_law_mask = T.
          !  These cells are now assigned to the calving domain.

          do j = 1, ny
             do i = 1, nx
                if (color(i,j) == fill_color) then
                   calving_domain_mask(i,j) = .true.
                else
                   calving_domain_mask(i,j) = .false.
                endif
             enddo
          enddo

          deallocate(color)

          ! Note: For this option, all cells in the calving domain have calving_law = T, so the logic below
          !  (calving_law_mask = T .and. calving_domain_mask = T) is redundant, but it does no harm.

       endif   ! calving_domain

       ! Calve ice where calving_law_mask = T and calving_domain_mask = T
       do j = 1, ny
          do i = 1, nx
             if (calving_law_mask(i,j) .and. calving_domain_mask(i,j)) then

                if (verbose_calving .and. this_rank==rtest .and. thck(i,j) > 0.0d0) then
                   print*, 'Calve ice: task, i, j, calving_thck =', this_rank, i, j, float_fraction_calve * thck(i,j)*thk0
                endif

                calving_thck(i,j) = float_fraction_calve * thck(i,j)
                thck(i,j) = thck(i,j) - calving_thck(i,j)
                !WHL TODO - Also handle tracers?  E.g., set damage(:,i,j) = 0.d0?
 
            endif
          enddo
       enddo

    endif   ! which_calving

    !WHL - debug
!    if (verbose_calving .and. main_task) then
!       j = jtest
!       print*, 'Calved ice: j =', j
!       print*, 'i, relx, topg, thck, usfc, calving_law_mask, calving_domain_mask:'
!       do i = iplot1, iplot2
!          print*, i, relx(i,j)*thk0, topg(i,j)*thk0, thck(i,j)*thk0, (topg(i,j) + thck(i,j))*thk0, &
!                  calving_law_mask(i,j), calving_domain_mask(i,j)
!       enddo
!    endif

    ! Remove any floating ice islands.
    ! Typically these will be removed by the calving scheme above, but if not, 
    !  then they need to be removed before calling the velocity solver,
    !  which will have problems in regions without any grounded ice.

    !WHL - debug - set up an ice island
!!    i = nx-10
!!    thck(i,:) = 0.0d0

    if (remove_floating_islands) then

       if (verbose_calving .and. main_task) then
          print*, 'Remove floating islands'
       endif

       ! Set the thickness thresholds for cells to be counted as active.
       ! Floating regions will be identified as ice islands unless they are connected to grounded ice
       !  along a path consisting only of active cells (i.e., cells with thk > thklim).
       ! Note: A limit of 0.0 does not work because it erroneously counts very thin floating cells as active.
       !       Then the algorithm can fail to identify floating regions that should be removed.
       ! Note: Optionally, we can specify a larger thickness threshold for floating ice, called
       !       floating_path_minthck. In other words, for a floating cell or region not to be removed,
       !       it must be connected to grounded ice via a path consisting of active grounded cells
       !       and/or active floating cells thicker than floating_path_minthck. 

       call glissade_remove_floating_islands(&
            thck,          relx,                  &
            topg,          eus,                   &
            thklim,        floating_path_minthck, &
            calving_thck)

    endif

    ! cleanup
    deallocate (ice_mask)
    deallocate (floating_mask)
    deallocate (ocean_mask)
    if (which_calving == CALVING_THCK_THRESHOLD) then
       deallocate (thick_or_grounded_mask)
       deallocate (protected_floating_mask)
    endif

    deallocate (calving_law_mask)
    deallocate (calving_domain_mask)

  end subroutine glissade_calve_ice

!---------------------------------------------------------------------------

  subroutine glissade_remove_floating_islands(&
       thck,          relx,                  &
       topg,          eus,                   &
       thklim,        floating_path_minthck, &
       calving_thck)

    ! Remove any floating ice islands. 
        
    ! The method is as follows: Initialize each cell to have either the initial color
    !  (if active ice is present) or the boundary color (if no ice is present).
    !  "Active" means that thck > thklim.
    ! Then loop through the cells on the processor. For each active grounded ice cell,
    !  assign the fill color and then recursively assign the fill color to any
    !  active cells with which it is connected (i.e., it shares an edge).
    ! Repeat the loop several times to allow communication between adjacent
    !  processors via halo updates.
    ! Any active cells that still have the initial color are floating ice islands.
    ! Remove this ice and add it to the calving field.

    use glissade_masks
    use glimmer_paramets, only: thk0

    real(dp), dimension(:,:), intent(inout) :: thck     !> ice thickness
    real(dp), dimension(:,:), intent(in)    :: relx     !> relaxed bedrock topography
    real(dp), dimension(:,:), intent(in)    :: topg     !> present bedrock topography
    real(dp), intent(in)      :: eus                    !> eustatic sea level
    real(dp), intent(in)      :: thklim                 !> minimum thickness for dynamically active ice
    real(dp), intent(in)      :: floating_path_minthck  !> minimum thickness along path from floating ice back to grounded ice
    real(dp), dimension(:,:), intent(inout) :: calving_thck   !> thickness lost due to calving in each grid cell;
                                                              !> on output, includes ice in floating islands

    integer :: nx, ny      ! horizontal grid dimensions

    integer :: i, j
    integer :: count, maxcount_fill  ! loop counters

    integer,  dimension(:,:), allocatable   ::  &
         ice_mask,          & ! = 1 where ice is present (thck > thklim), else = 0
         floating_mask,     & ! = 1 where ice is present and floating, else = 0
         color                ! integer 'color' for filling the calving domain (with CALVING_DOMAIN_OCEAN_CONNECT)

    !WHL - debug
    real(dp) :: sum_fill_local, sum_fill_global

    nx = size(thck,1)
    ny = size(thck,2)

    allocate (ice_mask(nx,ny))
    allocate (floating_mask(nx,ny))
    allocate (color(nx,ny))

    ! calculate masks
    ! Note: Passing in thklim = 0.0 does not work because it erroneously counts very thin floating cells as active.
    !       Then the algorithm can fail to identify floating regions that should be removed
    !       (since they are separated from any cells that are truly active).

     call glissade_get_masks(nx,            ny,            &
                             thck,          topg,          &
                             eus,           thklim,        &
                             ice_mask,      floating_mask)

    ! initialize
    ! Assign the initial color to cells with active ice and the boundary color to cells without active ice.

    do j = 1, ny
       do i = 1, nx
          if (ice_mask(i,j) == 1) then
             color(i,j) = initial_color
          else
             color(i,j) = boundary_color
          endif
       enddo
    enddo

    ! Loop through cells, identifying active cells with grounded ice.
    ! Fill each grounded cell and then recursively fill neighbor cells, whether grounded or not.
    ! We may have to do this several times to incorporate connections between neighboring processors.

    ! Note: Calling glissade_fill with the optional 'thck' and 'floating_path_minthck' arguments
    !       means that the cell must have thck > floating_path_minthck in order for neighbor cells 
    !       to be filled recursively.

    maxcount_fill = max(ewtasks,nstasks)

    if (verbose_calving .and. main_task) then
       print*, 'maxcount_fill =', maxcount_fill
    endif

    do count = 1, maxcount_fill

       if (count == 1) then   ! identify grounded cells that can seed the fill

          do j = 1, ny
             do i = 1, nx
                if (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) then  ! grounded ice
                   if (color(i,j) /= boundary_color .and. color(i,j) /= fill_color) then
                      ! assign the fill color to this cell, and recursively fill neighbor cells
                      call glissade_fill(nx, ny, i, j, color, thck, floating_path_minthck)
                   endif
                endif
             enddo
          enddo

       else  ! count > 1; first check for halo cells that were just filled on neighbor processors

          call parallel_halo(color)

          ! west halo layer
          i = nhalo
          do j = 1, ny
             if (color(i,j) == fill_color) call glissade_fill(nx, ny, i+1, j, color, thck, floating_path_minthck)
          enddo

          ! east halo layers
          i = nx - nhalo + 1
          do j = 1, ny
             if (color(i,j) == fill_color) call glissade_fill(nx, ny, i-1, j, color, thck, floating_path_minthck)
          enddo

          ! south halo layer
          j = nhalo
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color) call glissade_fill(nx, ny, i, j+1, color, thck, floating_path_minthck)
          enddo

          ! north halo layer
          j = ny-nhalo+1
          do i = nhalo+1, nx-nhalo  ! already checked halo corners above
             if (color(i,j) == fill_color) call glissade_fill(nx, ny, i, j-1, color, thck, floating_path_minthck)
          enddo

       endif  ! count = 1

       sum_fill_local = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             if (color(i,j) == fill_color) sum_fill_local = sum_fill_local + 1
          enddo
       enddo

       !WHL - If running a large problem, may want to reduce the frequency of this global sum
       sum_fill_global = parallel_reduce_sum(sum_fill_local)

       if (verbose_calving .and. main_task) then
          print*, 'this_rank, sum_fill_local, sum_fill_global:', this_rank, sum_fill_local, sum_fill_global
       endif

    enddo  ! count

    ! Any active cells that still have the initial color are part of floating ice islands.
    ! Remove ice in these cells, adding it to the calving field.
    do j = 1, ny
       do i = 1, nx
          if (color(i,j) == initial_color) then
             !WHL - debug
             if (verbose_calving .and. thck(i,j) > 0.0) &
                  print*, 'Remove floating island: task, i, j, thck =', this_rank, i, j, thck(i,j)*thk0
             calving_thck(i,j) = thck(i,j)
             thck(i,j) = 0.0d0
             !TODO - Also handle tracers?  E.g., set damage(:,i,j) = 0.d0?
          endif
       enddo
    enddo

    ! cleanup
    deallocate (ice_mask)
    deallocate (floating_mask)
    deallocate (color)

  end subroutine glissade_remove_floating_islands

!****************************************************************************

  recursive subroutine glissade_fill(nx,  ny,   &
                                     i,   j,    &
                                     color,     &
                                     thck,      &
                                     floating_path_minthck)

    ! Given a domain with an initial color, a boundary color and a fill color,
    ! assign the fill color to all cells that either (1) are prescribed to have
    ! the fill color or (2) are connected to cells with the fill color.

    integer, intent(in) :: nx, ny                       ! domain size
    integer, intent(in) :: i, j                         ! horizontal indices of current cell

    integer, dimension(nx,ny), intent(inout) :: &
         color                                          ! color (initial, fill or boundary)

    real(dp), dimension(nx,ny), intent(in), optional :: &
         thck                                           ! ice thickness

    real(dp), intent(in), optional :: &
         floating_path_minthck     ! min thickness along path from floating ice back to grounded ice

    !Note: If a cell is active and floating, but thck < floating_path_minthck,
    !       then fill this cell but do not call glissade_fill recursively.
    !      This is a method of removing long peninsulas of active but thin ice
    !       (thklim < thck < floating_path_minthck), which can have unrealistically
    !       large ice speeds.

    if (color(i,j) /= fill_color .and. color(i,j) /= boundary_color) then

       ! assign the fill color to this cell
       color(i,j) = fill_color

       if (present(thck) .and. present(floating_path_minthck)) then
          if (thck(i,j) < floating_path_minthck) then
             print*, 'Skip recursion: i, j, thck, minthck =', i, j, thck(i,j), floating_path_minthck
             return    ! skip the recursion
          endif
       endif

       ! recursively call this subroutine for each neighbor to see if it should be filled
       if (i > 1)  call glissade_fill(nx, ny, i-1, j,   color, thck, floating_path_minthck)
       if (i < nx) call glissade_fill(nx, ny, i+1, j,   color, thck, floating_path_minthck)
       if (j > 1)  call glissade_fill(nx, ny, i,   j-1, color, thck, floating_path_minthck)
       if (j < ny) call glissade_fill(nx, ny, i,   j+1, color, thck, floating_path_minthck)

    endif

  end subroutine glissade_fill

!---------------------------------------------------------------------------

end module glissade_calving

!---------------------------------------------------------------------------
