!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_basal_water.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glissade_basal_water

   use glimmer_global, only: dp
   use glimmer_paramets, only: eps11, eps08
   use glimmer_physcon, only: rhoi, rhow, grav, scyr
   use glimmer_log
   use glide_types
   use glide_diagnostics, only: point_diag
   use cism_parallel, only: main_task, this_rank, nhalo, parallel_type, parallel_halo

   implicit none

   private
   public :: glissade_basal_water_init, glissade_calcbwat, glissade_bwat_flux_routing

!!   logical, parameter :: verbose_bwat = .false.
   logical, parameter :: verbose_bwat = .true.

contains

!==============================================================

  subroutine glissade_basal_water_init(model)

    use glimmer_paramets, only: thk0

    type(glide_global_type) :: model

    select case (model%ho_options%which_ho_bwat)

    ! HO_BWAT_NONE:         basal water depth = 0
    ! HO_BWAT_CONSTANT:     basal water depth = prescribed constant
    ! HO_BWAT_LOCAL_TILL:   local basal till model with prescribed drainage rate
    ! HO_BWAT_FLUX_ROUTING: steady-state water routing with flux calculation

    case(HO_BWAT_CONSTANT)

       ! Set a constant water thickness where ice is present
       where (model%geometry%thck > model%numerics%thklim)
          model%basal_hydro%bwat(:,:) = model%basal_hydro%const_bwat / thk0
       elsewhere
          model%basal_hydro%bwat(:,:) = 0.0d0
       endwhere

    case default

       ! currently nothing to do for other options

    end select

  end subroutine glissade_basal_water_init

!==============================================================

  subroutine glissade_calcbwat(which_ho_bwat,    &
                               basal_hydro,      &
                               dt,               &
                               thck,             &
                               thklim,           &
                               bmlt,             &
                               bwat)

    ! Driver for updating basal water
    ! Note: This subroutine assumes SI units.
    ! Currently, only a few simple options are supported.

    use glide_types

    integer, intent(in) :: &
         which_ho_bwat     !> basal water options

    type(glide_basal_hydro), intent(inout) :: basal_hydro ! basal hydro object

    real(dp), intent(in) :: &
         dt,             & !> time step (s) 
         thklim            !> threshold for dynamically active ice (m)

    real(dp), dimension(:,:), intent(in) ::  &
         thck,           & !> ice thickness (m)
         bmlt              !> basal melt rate (m/s of ice)

    real(dp), dimension(:,:), intent(inout) ::  &
         bwat              !> basal water depth (m of water)

    ! local variables

    integer :: nx, ny    ! horizontal grid dimensions
    integer :: i, j

    real(dp) ::  &
         dbwat_dt        ! rate of change of bwat (m/s of water)

    select case (which_ho_bwat)

    ! HO_BWAT_NONE:         basal water depth = 0
    ! HO_BWAT_CONSTANT:     basal water depth = prescribed constant
    ! HO_BWAT_LOCAL_TILL:   local basal till model with prescribed drainage rate
    ! HO_BWAT_FLUX_ROUTING: steady-state water routing with flux calculation (handled in a different subroutine)

    case(HO_BWAT_NONE)

       bwat(:,:) = 0.0d0

    case(HO_BWAT_CONSTANT)

       ! Use a constant water thickness where ice is present, to force Tbed = Tpmp

       where (thck > thklim)
          bwat(:,:) = basal_hydro%const_bwat
       elsewhere
          bwat(:,:) = 0.0d0
       endwhere

     case(HO_BWAT_LOCAL_TILL)

        ! Change bwat based on input bmlt, while also draining water at a constant rate.
        ! Note: bmlt > 0 for ice melting, which increases bwat. Freeze-on (bmlt < 0) will reduce bwat.

        nx = size(bwat,1)
        ny = size(bwat,2)

        do j = 1, ny
           do i = 1, nx

              ! compute new bwat, given source (bmlt) and sink (drainage)
              dbwat_dt = bmlt(i,j)*rhoi/rhow - basal_hydro%c_drainage/scyr  ! convert c_drainage from m/yr to m/s
              bwat(i,j) = bwat(i,j) + dbwat_dt*dt

              ! limit to the range [0, bwat_till_max]
              bwat(i,j) = min(bwat(i,j), basal_hydro%bwat_till_max)
              bwat(i,j) = max(bwat(i,j), 0.0d0)

           enddo
        enddo

    end select

  end subroutine glissade_calcbwat

!==============================================================

  subroutine glissade_bwat_flux_routing(&
       nx,            ny,            &
       dx,            dy,            &
       parallel,                     &
       itest, jtest,  rtest,         &
       flux_routing_scheme,          &
       thck,          topg,          &
       thklim,                       &
       bwat_mask,     floating_mask, &
       bmlt_hydro,                   &
       bwatflx,       bwat_diag,     &
       head,          grad_head)

    ! Compute the subglacial water flux and water depth using a steady-state flux routing scheme.
    ! Water is routed down the hydropotential. For routing purposes, assume p_w = p_i (i.e., N = 0).
    ! Elsewhere in the code, N can be computed as a function of the water depth diagnosed here.

    !TODO - Pass in a potential for basal freezing (bfrz_pot)?
    !       Return the actual bfrz.

    use cism_parallel, only: tasks   ! while code is serial only

    ! Input/output arguments

    integer, intent(in) ::  &
         nx, ny,             &      ! number of grid cells in each direction
         itest, jtest, rtest        ! coordinates of diagnostic point

    real(dp), intent(in) ::  &
         dx, dy                     ! grid cell size (m)

    type(parallel_type), intent(in) ::  &
         parallel                   ! info for parallel communication

    integer, intent(in) ::  &
         flux_routing_scheme        ! flux routing scheme: D8, Dinf or FD8; see subroutine route_basal_water

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                    & ! ice thickness (m)
         topg,                    & ! bed topography (m)
         bmlt_hydro                 ! meltwater input rate (m/s); can include an englacial or surface source

    real(dp), intent(in) ::  &
         thklim                     ! minimum ice thickness for basal melt and hydropotential calculations (m)
                                    ! Note: This is typically model%geometry%thklim_temp

    integer, dimension(nx,ny), intent(in) ::  &
         bwat_mask,               & ! mask to identify cells through which basal water is routed;
                                    ! = 0 for floating and ocean cells; cells at global domain edge;
                                    !  and cells with thck = 0 and forced negative SMB
         floating_mask              ! = 1 if ice is present (thck > thklim) and floating, else = 0

    ! Note: The field bwat_diag is the steady-state basal water depth, diagnosed as a function of bwatflx.
    !       Outside this subroutine, bwat_diag can be used to compute effective pressure,
    !        but it does not carry any enthalpy, so it is ignored in glissade_therm.
    !       The bwat field, on the other hand, is prognosed in the local basal till model;
    !        this water must be frozen in glissade_therm before the bed can cool below Tpmp.
    ! Note: The hydropotential phi = rhow*grav*head.
    !       We follow Sommers et al. (2018) in using 'head', which has convenient units of meters.

    real(dp), dimension(nx,ny), intent(out) ::  &
         bwatflx,                 & ! basal water flux (m/s)
         bwat_diag,               & ! diagnosed basal water depth (m)
         head,                    & ! hydraulic head (m)
         grad_head                  ! gradient of hydraulic head (m/m), averaged to cell centers

    ! Local variables

    integer :: i, j, p

    real(dp), dimension(nx, ny) ::  &
         lakes                      ! difference between filled head and original head (m)

    ! parameters related to subglacial fluxes
    ! The water flux q is given by Sommers et al. (2018), Eq. 5:
    !
    !           q = (b^3*g)/[(12*nu)*(1 + omega*Re)] * (-grad(h))
    !
    ! where q = basal water flux per unit width (m^2/s) = bwatflx/dx
    !       b = water depth (m) = bwat
    !       g = gravitational constant (m/s^2) = grad
    !      nu = kinematic viscosity of water (m^2/s) = visc_water
    !   omega = parameter controlling transition between laminar and turbulent flow
    !      Re = Reynolds number (large for turbulent flow)
    !       h = hydraulic head (m)
    !
    ! Note: In the equation above and the calculation below, bwatflx has units of m^3/s,
    !        i.e., volume per second entering and exiting a grid cell.
    !       For output, bwatflx has units of m/s, i.e. volume per unit area per year entering and exiting a grid cell.
    !       With the latter convention, bwatflx is independent of grid resolution.
    !
    ! By default, we set Re = 0, which means the flow is purely laminar, as in Sommers et al. (2018), Eq. 6.
    !TODO - Compute Re based on the flux.

    real(dp), parameter ::  &
         visc_water = 1.787e-6,            & ! kinematic viscosity of water (m^2/s); Sommers et al. (2018), Table 2
         omega_hydro = 1.0d-3                ! omega (unitless) in Sommers et al (2018), Eq. 6

    real(dp), parameter ::  &
         p_flux_to_depth = 2.0d0,          & ! exponent for water depth; = 2 if q is proportional to b^3
         q_flux_to_depth = 1.0d0             ! exponent for potential gradient; = 1 if q is linearly proportional to grad(h)

    real(dp) :: c_flux_to_depth              ! proportionality coefficient in Sommers et al., Eq. 6
    real(dp) :: Reynolds                     ! Reynolds number (unitless), = 0 for pure laminar flow

    integer :: nx_test, ny_test
    real(dp), dimension(:,:), allocatable :: phi_test
    integer,  dimension(:,:), allocatable :: mask_test

    if (verbose_bwat .and. this_rank == rtest) then
       write(6,*) 'In glissade_bwat_flux_routing: rtest, itest, jtest =', rtest, itest, jtest
    endif

    ! Uncomment if the following fields are not already up to date in halo cells
!    call parallel_halo(thk,  parallel)
!    call parallel_halo(topg, parallel)
    call parallel_halo(bmlt_hydro, parallel)
    !TODO - Add bfrz?


    ! Compute the hydraulic head
    ! For purposes of flux routing, assume N = 0.
    ! Then the head depends only on basal topography and ice thickness.

    call compute_head(&
         nx,     ny,    &
         thck,          &
         topg,          &
         thklim,        &
         floating_mask, &
         head)

    if (verbose_bwat) then
       call point_diag(thck,      'thck (m)',    itest, jtest, rtest, 7, 7)
       call point_diag(topg,      'topg (m)',    itest, jtest, rtest, 7, 7)
       call point_diag(bmlt_hydro*scyr, 'bmlt_hydro (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(bwat_mask, 'bwat_mask',   itest, jtest, rtest, 7, 7)
       call point_diag(head,      'Before fill, head (m)', itest, jtest, rtest, 7, 7)
    endif

    ! Route basal water down the gradient of hydraulic head, giving a water flux
    ! TODO - Pass in bfrz_pot, return bfrz?

    call route_basal_water(&
         nx,      ny,            &
         dx,      dy,            &
         parallel,               &
         itest, jtest, rtest,    &
         flux_routing_scheme,    &
         head,                   &
         bmlt_hydro,             &
         bwat_mask,              &
         bwatflx,                &
         lakes)

    call parallel_halo(bwatflx, parallel)

    ! Convert the water flux to a basal water depth
    ! Note: bwat_diag is treated differently from bwat, for the following reason.
    ! In the thermal solve, the basal temperature is held at Tpmp wherever bwat > 0.
    ! This is appropriate when bwat is prognosed from local basal melting.
    ! For the flux-routing scheme, however, we can diagnose nonzero bwat beneath ice
    !  that is frozen to the bed with T < Tpmp (due to basal melting upstream).
    ! If passed to the thermal solver, this bwat could drive a sudden large increase in basal temperature.

    ! Set parameters in the flux-to-depth relationship
    !TODO: Set Reynolds as a function of qflx and visc_water

    Reynolds = 0.0d0
    c_flux_to_depth = 1.0d0/((12.0d0*visc_water)*(1.0d0 + omega_hydro*Reynolds))

    call flux_to_depth(&
         nx,       ny,           &
         dx,       dy,           &
         itest,    jtest, rtest, &
         bwatflx,                &
         head,                   &
         c_flux_to_depth,        &
         p_flux_to_depth,        &
         q_flux_to_depth,        &
         bwat_mask,              &
         bwat_diag,              &
         grad_head)

    call parallel_halo(bwat_diag, parallel)
    call parallel_halo(grad_head, parallel)

    ! Convert bwatflx units from m^3/s to m/s for output
    bwatflx(:,:) = bwatflx(:,:) / (dx*dy)

    if (verbose_bwat) then
       call point_diag(bwatflx*scyr, 'Final bwatflx (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(bwat_diag, 'Diagnosed bwat (m)',  itest, jtest, rtest, 7, 7, '(f10.6)')
    endif

  end subroutine glissade_bwat_flux_routing

!==============================================================

  subroutine compute_head(&
       nx,      ny,   &
       thck,          &
       topg,          &
       thklim,        &
       floating_mask, &
       head)

    !  Approximate the hydraulic head as the bed elevation plus the scaled water pressure:
    !
    !     head = z_b + p_w / (rhow*g)
    !
    !  where z_b = bed elevation (m) = topg
    !        p_w = water pressure (Pa) = p_i - N
    !        p_i = ice overburden pressure = rhoi*g*H
    !          N = effective pressure (Pa) = part of overburden not supported by water
    !          H = ice thickness (m)
    !
    !  If we make the approximation p_w = p_i = rhoi*g*H (i.e., N = 0), then
    !
    !     head = z_b + (rhoi/rhow) * H
    !
    !  Note: Elsewhere, in glissade_calc_effecpress, N is computed explicitly and is not equal to 0.

    implicit none

    ! Input/output variables

    integer, intent(in) ::  &
         nx, ny                   ! number of grid cells in each direction

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                  & ! ice thickness (m)
         topg                     ! bed elevation (m)

    real(dp), intent(in) ::  &
         thklim                   ! minimum ice thickness for bmlt and head calculations

    integer, dimension(nx,ny), intent(in) ::  &
         floating_mask            ! = 1 if ice is present (thck > thklim) and floating, else = 0

    real(dp), dimension(nx,ny), intent(out) ::  &
         head                     ! hydraulic head (m)

    where (thck > thklim .and. floating_mask /= 1)
       head = topg + (rhoi/rhow)*thck
    elsewhere
       head = max(topg, 0.0d0)
    endwhere

  end subroutine compute_head

!==============================================================

  subroutine route_basal_water(&
         nx,         ny,         &
         dx,         dy,         &
         parallel,               &
         itest, jtest, rtest,    &
         flux_routing_scheme,    &
         head,                   &
         bmlt_hydro,             &
         bwat_mask,              &
         bwatflx,                &
         lakes)

    ! Route water from the basal melt field to its destination, recording the water flux along the way.
    ! Water flow direction is determined according to the gradient of the hydraulic head.
    ! For the algorithm to work correctly, surface depressions must be filled,
    !  so that all cells have an outlet to the ice sheet margin.
    ! This results in the lakes field, which is the difference between the filled head and the original head.
    !  The method used is by Quinn et. al. (1991).
    !
    ! Based on code by Jesse Johnson (2005), adapted from the glimmer_routing file by Ian Rutt.

    ! TODO - Pass in bfrz_pot, return bfrz.

    use cism_parallel, only: parallel_global_sum

    !WHL - debug
    use cism_parallel, only: parallel_globalindex, parallel_reduce_max

    implicit none

    integer, intent(in) ::  &
         nx, ny,               & ! number of grid cells in each direction
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), intent(in) ::  &
         dx, dy                  ! grid cell size (m)

    type(parallel_type), intent(in) ::  &
         parallel                ! info for parallel communication

    integer, intent(in) ::  &
         flux_routing_scheme     ! flux routing scheme: D8, Dinf or FD8

    real(dp), dimension(nx,ny), intent(in)  ::  &
         bmlt_hydro              ! meltwater input rate (m/s)

    real(dp), dimension(nx,ny), intent(inout)  ::  &
         head                    ! hydraulic head (m)
                                 ! intent inout because it can be filled and adjusted below

    integer, dimension(nx,ny), intent(in)  ::  &
         bwat_mask               ! mask to identify cells through which basal water is routed;
                                 ! = 1 where ice is present and not floating

    real(dp), dimension(nx,ny), intent(out) ::  &
         bwatflx,             &  ! water flux through a grid cell (m^3/s)
         lakes                   ! lakes field, difference between filled and original head

    ! Local variables

    integer :: nlocal            ! number of locally owned cells
    integer :: count, count_max  ! iteration counters
    integer :: i, j, k, ii, jj, ip, jp, p
    integer :: i1, j1, i2, j2, itmp, jtmp

    logical :: finished    ! true when an iterative loop has finished

    integer,  dimension(:,:), allocatable ::  &
         sorted_ij         ! i and j indices of all cells, sorted from low to high values of head

    real(dp), dimension(-1:1,-1:1,nx,ny) ::  &
         flux_fraction, &  ! fraction of flux from each cell that flows downhill to each of 8 neighbors
         bwatflx_halo      ! water flux (m^3/s) routed to a neighboring halo cell; routed further in next iteration

    real(dp), dimension(nx,ny) ::  &
         head_filled,   &  ! head after depressions are filled (m)
         bwatflx_accum, &  ! water flux (m^3/s) accumulated over multiple iterations
         sum_bwatflx_halo  ! bwatflx summed over the first 2 dimensions in each grid cell

    integer, dimension(nx,ny) ::  &
         local_mask,     & ! = 1 for cells owned by the local processor, else = 0
         halo_mask,      & ! = 1 for the layer of halo cells adjacent to locally owned cells, else = 0
         margin_mask       ! = 1 for cells at the grounded ice margin, as defined by bwat_mask, else = 0

    real(dp) :: &
         total_flux_in,  & ! total input flux (m^3/s), computed as sum of bmlt_hydro*dx*dy
         total_flux_out, & ! total output flux (m^3/s), computed as sum of bwatflx at ice margin
         err,            & ! water conservation error
         global_flux_sum   ! flux sum over all cells in global domain

    character(len=100) :: message

    !WHL - debug
    real(dp) :: bmlt_max, bmlt_max_global
    integer :: imax, jmax, rmax, iglobal, jglobal
    ! Allocate the sorted_ij array

    nlocal = parallel%own_ewn * parallel%own_nsn
    allocate(sorted_ij(nlocal,2))

    ! Compute mask of locally owned and halo cells.
    ! These masks are used to transfer fluxes between processors on subsequent iterations.

    local_mask = 0
    halo_mask = 0
    do j = nhalo, ny-nhalo+1
       do i = nhalo, nx-nhalo+1
          if (j == nhalo .or. j == ny-nhalo+1 .or. i == nhalo .or. i == nx-nhalo+1) then
             halo_mask(i,j) = 1
          elseif (j > nhalo .or. j <= ny-nhalo .or. i > nhalo .or. i <= nx-nhalo+1) then
             local_mask(i,j) = 1
          endif
       enddo
    enddo

    ! Fill depressions in the head field, so that no interior cells are sinks

    call fill_depressions(&
         nx,    ny,            &
         parallel,             &
         itest, jtest, rtest,  &
         head,                 &
         bwat_mask,            &
         head_filled)

    ! Compute the lake depth
    lakes = head_filled - head

    ! Update head with the filled values
    head = head_filled

    if (verbose_bwat) then
       call point_diag(head,  'head (m)',  itest, jtest, rtest, 7, 7)
       call point_diag(lakes, 'lakes (m)', itest, jtest, rtest, 7, 7)
    endif

    ! Sort heights.
    ! The sorted_ij array stores the i and j index for each locally owned cell, from lowest to highest value.

    call sort_heights(&
         nx,    ny,    nlocal,  &
         itest, jtest, rtest,   &
         head,  sorted_ij)

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'sorted, from the top:'
       do k = nlocal, nlocal-10, -1
          i = sorted_ij(k,1)
          j = sorted_ij(k,2)
          print*, k, i, j, head(i,j)
       enddo
    endif

    call get_flux_fraction(&
         nx,    ny,    nlocal,  &
         dx,    dy,             &
         itest, jtest, rtest,   &
         flux_routing_scheme,   &
         sorted_ij,             &
         head,                  &
         bwat_mask,             &
         flux_fraction)

    ! Initialize bwatflx in locally owned cells with the basal melt, which will be routed downslope.
    ! Multiply by area, so units are m^3/s.
    ! The halo water flux, bwatflx_halo, holds water routed to halo cells;
    !  it will be routed downhill on the next iteration.
    ! The accumulated flux, bwatflx_accum, holds the total flux over multiple iterations.
    ! Note: This subroutine conserves water only if bmlt_hydro >= 0 everywhere.
    !       One way to account for refreezing would be to do the thermal calculation after
    !        computing bwat in this subroutine.  At that point, refreezing would take away
    !        from the bwat computed here.  In the next time step, positive values of bmlt_hydro
    !        would provide a new source for bwat.
    ! In other words, the sequence would be:
    ! (1) Ice transport and calving
    ! (2) Basal water routing: apply bmlt_hydro and diagnose bwat
    ! (3) Vertical heat flow:
    !     (a) compute bmlt_hydro
    !     (b) use bmlt_hydro < 0 to reduce bwat
    !     (c) save bmlt_hydro > 0 for the next time step (and write to restart)
    ! (4) Diagnose velocity

    bwatflx = 0.0d0
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          bwatflx(i,j) = bmlt_hydro(i,j) * dx * dy
          bwatflx(i,j) = max(bwatflx(i,j), 0.0d0)   ! not conservative unless refreezing is handled elsewhere
       enddo
    enddo
    bwatflx_halo = 0.0d0
    bwatflx_accum = 0.0d0

    ! Compute total input of meltwater (m^3/s)
    total_flux_in = parallel_global_sum(bwatflx, parallel)

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'Total input basal melt flux (m^3/s):', total_flux_in
    endif

    ! Loop over locally owned cells, from highest to lowest.
    ! During each iteration, there are two possible outcomes for routing:
    ! (1) Routed to the ice sheet margin, to a cell with bwat_mask = 0.
    !     In this case, the routing of that flux is done.
    ! (2) Routed to a halo cell, i.e. a downslope cell on a neighboring processor.
    !     In this case, the flux will be routed further downhill on the next iteration.
    ! When all the water has been routed to the margin, we are done.

    count = 0
    ! Note: It is hard to predict how many iterations will be sufficient.
    !       With Dinf or FD8, we can have flow back and forth across processor boundaries,
    !        requiring many iterations to reach the margin.
    !       For Greenland 4 km, Dinf requires ~20 iterations on 4 cores, and FD8 can require > 40.
    !       For Antarctica 8 km, FD8 can require > 50.
    count_max = 100
    finished = .false.

    do while (.not.finished)

       count = count + 1

       if (verbose_bwat .and. this_rank == rtest) then
          print*, 'flux routing, count =', count
       endif

       do k = nlocal, 1, -1

          ! Get i and j indices of current point
          i = sorted_ij(k,1)
          j = sorted_ij(k,2)

          ! Apportion the flux among downslope neighbors
          if (bwat_mask(i,j) == 1 .and. bwatflx(i,j) > 0.0d0) then
             do jj = -1,1
                do ii = -1,1
                   ip = i + ii
                   jp = j + jj
                   if (flux_fraction(ii,jj,i,j) > 0.0d0) then
                      if (halo_mask(ip,jp) == 1) then
                         bwatflx_halo(ii,jj,i,j) = bwatflx(i,j)*flux_fraction(ii,jj,i,j)
                         if (verbose_bwat .and. this_rank==rtest .and. i==itest .and. j==jtest .and. count <= 2) then
                            print*, 'Flux to halo, i, j, ii, jj, flux:', &
                                 i, j, ii, jj, bwatflx(i,j)*flux_fraction(ii,jj,i,j)
                         endif
                      elseif (local_mask(ip,jp) == 1) then
                         bwatflx(ip,jp) = bwatflx(ip,jp) + bwatflx(i,j)*flux_fraction(ii,jj,i,j)
                      endif
                   endif   ! flux_fraction > 0
                enddo
             enddo
          endif

       enddo  ! loop from high to low

       ! Accumulate bwatflx from the latest iteration.
       ! Reset to zero for the next iteration, if needed.

       bwatflx_accum = bwatflx_accum + bwatflx
       bwatflx = 0.0d0

       if (verbose_bwat .and. this_rank == rtest .and. count <= 2) then
          i = itest
          j = jtest
          print*, 'i, j, bwatflx_accum:', i, j, bwatflx_accum(i,j)
       endif

       ! If bwatflx_halo = 0 everywhere, then we are done.
       ! (If the remaining flux is very small (< eps11), discard it to avoid
       !  unnecessary extra iterations.)
       ! If bwatflx_halo remains, then communicate it to neighboring tasks and
       !  continue routing on the next iteration.

       do j = 1, ny
          do i = 1, nx
             sum_bwatflx_halo(i,j) = sum(bwatflx_halo(:,:,i,j))
!             if (verbose_bwat .and. sum_bwatflx_halo(i,j) > eps11 .and. count > 50) then
!               print*, 'Nonzero bwatflx_halo, count, rank, i, j, sum_bwatflx_halo:', &
!                     count, this_rank, i, j, sum_bwatflx_halo(i,j)
!               call parallel_globalindex(i, j, iglobal, jglobal, parallel)
!               print*, '     iglobal, jglobal:', iglobal, jglobal
!             endif
          enddo
       enddo
       global_flux_sum = parallel_global_sum(sum_bwatflx_halo, parallel)

       if (verbose_bwat .and. count <= 2) then
          if (this_rank == rtest) then
             write(6,*) 'Before halo update, sum of bwatflx_halo:', global_flux_sum
          endif
          call point_diag(sum_bwatflx_halo, 'sum_bwatflx_halo', itest, jtest, rtest, 7, 7)
       endif

       if (global_flux_sum > eps11) then

          finished = .false.

          ! Communicate bmltflx_halo to the halo cells of neighboring processors
          call parallel_halo(bwatflx_halo(:,:,:,:), parallel)

          ! bmltflx_halo is now available in the halo cells of the local processor.
          ! Route downslope to the adjacent locally owned cells.
          ! These fluxes will be routed further downslope during the next iteration.

          do j = 2, ny-1
             do i = 2, nx-1
                if (halo_mask(i,j) == 1 .and. sum(bwatflx_halo(:,:,i,j)) > 0.0d0) then
                   do jj = -1,1
                      do ii = -1,1
                         if (bwatflx_halo(ii,jj,i,j) > 0.0d0) then
                            ip = i + ii
                            jp = j + jj
                            if (local_mask(ip,jp) == 1) then
                               bwatflx(ip,jp) = bwatflx(ip,jp) + bwatflx_halo(ii,jj,i,j)
                               if (verbose_bwat .and. ip==itest .and. jp==jtest .and. this_rank==rtest &
                                    .and. count <= 2) then
                                  print*, 'Nonzero bwatflx from halo, rank, i, j:', &
                                       this_rank, ip, jp, bwatflx_halo(ii,jj,i,j)
                               endif
                            endif
                         endif   !  bwatflx_halo > 0 to a local cell
                      enddo   ! ii
                   enddo   ! jj
                endif   ! bwatflx_halo > 0 from this halo cell
             enddo   ! i
          enddo   ! j

          ! Reset bwatflx_halo for the next iteration
          bwatflx_halo = 0.0d0

          global_flux_sum = parallel_global_sum(bwatflx, parallel)
          if (verbose_bwat .and. this_rank == rtest .and. count <= 2) then
             ! Should be equal to the global sum of bwatflx_halo computed above
             print*, 'After halo update, sum(bwatflx from halo) =', global_flux_sum
             print*, ' '
          endif

       else   ! bwatflx_halo = 0 everywhere; no fluxes to route to adjacent processors
          if (verbose_bwat .and. this_rank == rtest) print*, 'Done routing fluxes'
          finished = .true.
          bwatflx = bwatflx_accum
       endif

       if (count > count_max) then
          call write_log('Hydrology error: too many iterations in route_basal_water', GM_FATAL)
       endif

    enddo  ! finished routing

    ! Identify cells just beyond the ice sheet margin, which can receive from upstream but not send downstream
    where (bwat_mask == 0 .and. bwatflx > 0.0d0)
       margin_mask = 1
    elsewhere
       margin_mask = 0
    endwhere

    ! Compute total output of meltwater (m^3/s) and check that input = output, within roundoff.

    total_flux_out = parallel_global_sum(bwatflx*margin_mask, parallel)

    if (verbose_bwat .and. this_rank == rtest) then
       print*, 'Total output basal melt flux (m^3/s):', total_flux_out
       print*, 'Difference between input and output =', total_flux_in - total_flux_out
    endif

    ! Not sure if a threshold of eps11 is large enough.  Increase if needed.
    if (total_flux_in > 0.0d0) then
       err = abs(total_flux_in - total_flux_out)
       if (err > eps11) then
!          write(message,*) 'Hydrology error: total water not conserved, error =', err
!          call write_log(message, GM_FATAL)
          write(message,*) 'WARNING: Hydrology error: total water not conserved, error =', err
          call write_log(message, GM_WARNING)
       endif
    endif

    ! clean up
    deallocate(sorted_ij)

  end subroutine route_basal_water

!==============================================================

  subroutine flux_to_depth(&
         nx,       ny,           &
         dx,       dy,           &
         itest,    jtest, rtest, &
         bwatflx,                &
         head,                   &
         c_flux_to_depth,        &
         p_flux_to_depth,        &
         q_flux_to_depth,        &
         bwat_mask,              &
         bwat_diag,              &
         grad_head)

    !  Assuming that the flow is steady state, this function simply solves
    !               flux = depth * velocity
    !  for the depth, assuming that the velocity is a function of depth,
    !  and pressure potential. This amounts to assuming a Weertman film,
    !  or Manning flow, both of which take the form of a constant times water
    !  depth to a power, times grad(head) to a power.
    !
    !  The subroutine also returns grad(head), which enters the flux calculation

    use glissade_grid_operators, only: glissade_gradient_at_edges

    ! Input/ouput variables

    integer, intent(in) ::  &
         nx, ny,               & ! number of grid cells in each direction
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), intent(in) ::  &
         dx, dy                   ! grid spacing in each direction

    real(dp), dimension(nx,ny), intent(in) :: &
         bwatflx,               & ! basal water flux (m^3/s)
         head                     ! hydraulic head (m)

    real(dp), intent(in) ::  &
         c_flux_to_depth,       & ! constant of proportionality
         p_flux_to_depth,       & ! exponent for water depth
         q_flux_to_depth          ! exponent for potential_gradient

    integer, dimension(nx,ny), intent(in)  ::  &
         bwat_mask                ! mask to identify cells through which basal water is routed;
                                  ! = 1 where ice is present and not floating

    real(dp), dimension(nx,ny), intent(out)::  &
         bwat_diag,             & ! water depth diagnosed from bwatflx
         grad_head                ! gradient of hydraulic head (m/m), averaged to cell centers

    ! Local variables

    integer :: i, j, p

    real(dp), dimension(nx-1,ny) ::  &
         dhead_dx                 ! gradient component on E edges

    real(dp), dimension(nx,ny-1) ::  &
         dhead_dy                 ! gradient component on N edges

    real(dp) ::  &
         dhead_dx_ctr, dhead_dy_ctr, & ! gradient components averaged to cell centers
         p_exponent                    ! p-dependent exponent in bwat expression

    integer, dimension(nx,ny) ::  &
         ice_mask                 ! mask passed to glissade_gradient_at edges; = 1 everywhere

    ice_mask = 1

    ! Compute gradient components at cell edges
    ! HO_GRADIENT_MARGIN_LAND: Use all field values when computing the gradient, including values in ice-free cells.

    call glissade_gradient_at_edges(&
         nx,       ny,       &
         dx,       dy,       &
         head,               &
         dhead_dx, dhead_dy, &
         ice_mask,           &
         gradient_margin_in = HO_GRADIENT_MARGIN_LAND)

    grad_head = 0.0d0  ! will remain 0 in outer row of halo cells
    do j = 2, ny-1
       do i = 2, nx-1
          dhead_dx_ctr = 0.5d0 * (dhead_dx(i-1,j) + dhead_dx(i,j))
          dhead_dy_ctr = 0.5d0 * (dhead_dy(i,j-1) + dhead_dy(i,j))
          grad_head(i,j) = sqrt(dhead_dx_ctr**2 + dhead_dy_ctr**2)
       enddo
    enddo

    !TODO - If a halo update is needed for grad_head, then pass in 'parallel'.  But may not be needed.
!!    call parallel_halo(grad_head, parallel)

    !WHL - debug
    if (verbose_bwat) then
       call point_diag(grad_head, 'grad_head (m/m)', itest, jtest, rtest, 7, 7)
    endif

    p_exponent = 1.d0 / (p_flux_to_depth + 1.d0)

    ! Note: In Sommers et al. (2018), Eq. 5, the basal water flux q (m^2/s) is
    !              q = (b^3 * g) / [(12*nu)(1 + omega*Re)] * (-grad(h))
    !       where nu = kinematic viscosity of water = 1.787d-06 m^2/s
    !          omega = 0.001
    !             Re = Reynolds number
    !
    ! Following Aleah's formulation:
    !          F = b^3 * c * g * dx * -grad(h) where c = 1/[(12*nu)(1 + omega*Re)]
    !        b^3 = F / [c * g * dx * -grad(h)]
    !          b = { F / [c * g * dx * -grad(h)] }^(1/3)
    !
    ! In the context of a formulation with general exponents,
    ! we have q_flux_to_depth = 1 and p_flux_to_depth = 2 (so p_exponent = 1/3)
    !
    ! Jesse's Glimmer code had this:
    !       bwat_diag = ( bwatflx / (c_flux_to_depth * scyr * dy * grad_wphi**q_flux_to_depth) ) ** exponent
    ! which is missing the grav term and seems to have an extra scyr term.
    ! Also, c_flux_to_depth = 1 / (12 * 1.6d-3) in Jesse's code.  Note exponent of d-3 instead of d-6 for nu.
    !
    ! Note: Assuming dx = dy
    ! TODO: Modify for the case dx /= dy?
    ! TODO: Compute c_flux_to_depth with a nonzero omega and Reynolds number

    ! Rescale this equation if using grad_phi_hydro (insert factor of rhow)
    where (grad_head /= 0.d0 .and. bwat_mask == 1)
       bwat_diag = ( bwatflx / (c_flux_to_depth * grav * dy * grad_head**q_flux_to_depth) ) ** p_exponent
    elsewhere
       bwat_diag = 0.d0
    endwhere

  end subroutine flux_to_depth

!==============================================================

  subroutine fill_depressions(&
       nx,    ny,            &
       parallel,             &
       itest, jtest, rtest,  &
       phi_in,               &
       phi_mask,             &
       phi)

    ! Fill depressions in the input field, phi_in.
    ! The requirements for the output field, phi_out, are:
    ! (1) phi_out >= phi_in everywhere
    ! (2) For each cell with phi_mask = 1, there is a descending path to the boundary.
    !     That is, phi1 >= phi2 for any two adjacent cells along the path, where the flow
    !     is from cell 1 to cell 2.
    ! (3) phi_out is the lowest surface consistent with properties (1) and (2).
    !
    ! The algorithm is based on this paper:
    ! Planchon, O., and F. Darboux (2001): A fast, simple and versatile algorithm
    !  to fill the depressions of digital elevation models, Catena (46), 159-176.
    !
    ! The basic idea is:
    ! Let phi = the current best guess for phi_out.
    ! Initially, set phi = phi_in on the boundary, and set phi = a large number elsewhere.
    ! Loop through the domain.  For each cell c, with value phi(c) not yet fixed as a known value,
    !  compute phi_min8(n), the current minimum of phi in the 8 neighbor cells.
    ! If phi_in(c) > phi_min8(n) + eps, then set phi(c) = phi_in(c) and mark that cell as having a known value,
    !  since phi(c) cannot go any lower.  Here, eps is a small number greater than zero.
    ! If phi_in(c) < phi_min8(n) + eps, but phi(c) > phi_min8(c) + eps, set phi(c) = phi_min8(n) + eps.
    !  Do not mark the cell as having a known value, because it might be lowered further.
    ! Continue until no further lowering of phi is possible.  At that point, phi = phi_out.
    ! Note: Setting eps = 0 would result in flat surfaces that would need to be fixed later.

    use cism_parallel, only: parallel_reduce_sum
    use cism_parallel, only: parallel_globalindex

    implicit none

    ! Input/output variables

    integer, intent(in) ::  &
         nx, ny,               & ! number of grid cells in each direction
         itest, jtest, rtest     ! coordinates of diagnostic point

    type(parallel_type), intent(in) ::  &
         parallel            ! info for parallel communication

    real(dp), dimension(nx,ny), intent(in) :: &
         phi_in              ! input field with depressions to be filled

    integer, dimension(nx,ny), intent(in) ::  &
         phi_mask            ! = 1 in the domain where depressions need to be filled, else = 0

    real(dp), dimension(nx,ny), intent(out) :: &
         phi                 ! output field with depressions filled

    ! Local variables --------------------------------------

    logical, dimension(nx,ny) ::  &
         known               ! = true for cells where the final phi(i,j) is known

    integer :: &
         local_lowered,    & ! local sum of cells where phi is lowered
         global_lowered      ! global sum of cells where phi is lowered

    real(dp) :: &
         phi_min8            ! current minval of phi in a cell's 8 neighbors,
                             ! considering only cells with phi_mask = 1

    real(dp) :: epsilon      ! small increment in phi, either epsilon_edge or epsilon_diag

    logical :: finished      ! true when an iterative loop has finished

    integer :: count         ! iteration counter

    integer :: i, j, ii, jj, ip, jp, p
    integer :: iglobal, jglobal
    integer :: i1, i2, istep, j1, j2, jstep

    real(dp), parameter :: big_number = 1.d+20   ! initial large value for phi

    ! According to Planchon & Darboux (2001), there should be one value of epsilon for edge neighbors
    !  and another value for corner neighbors.
    real(dp), parameter :: &
         epsilon_edge = 1.d-3,          & ! small increment in phi to avoid flat regions, applied to edge neighbors
         epsilon_diag = 1.d-3*sqrt(2.d0)  ! small increment in phi to avoid flat regions, applied to diagonal neighbors

    !WHL - Typically, it takes ~10 iterations to fill all depressions on a large domain.
    integer, parameter :: count_max = 100

!!    logical, parameter :: verbose_depression = .false.
    logical, parameter :: verbose_depression = .true.

    ! Initial halo updates, in case phi_in and phi_mask are not up to date in halo cells
    call parallel_halo(phi_in, parallel)
    call parallel_halo(phi_mask, parallel)

    ! Initialize phi to a large value
    where (phi_mask == 1)
       phi = big_number
       known = .false.
    elsewhere
       phi = 0.0d0
       known = .true.
    endwhere

    ! Set phi = phi_in for boundary cells (cells with phi_mask = 0).
    where (phi_mask == 0)
       phi = phi_in
       known = .true.
    endwhere

    if (verbose_depression) then
       call point_diag(phi, 'Initial phi', itest, jtest, rtest, 7, 7, '(es10.3)')
    endif

    count = 0
    finished = .false.

    do while (.not.finished)

       count = count + 1
       local_lowered = 0

       if (verbose_depression .and. this_rank == rtest) then
          write(6,*) ' '
          write(6,*) 'fill_depressions, count =', count
       endif

       ! Loop through cells
       ! Iterate until phi cannot be lowered further.
       !
       ! To vary the route through the cells and reduce the required number of iterations,
       ! we alternate between four possible sequences:
       ! (1) j lo to hi, i lo to hi
       ! (2) j hi to lo, i hi to lo
       ! (3) j lo to hi, i hi to lo
       ! (4) j hi to lo, i lo to hi
       ! Other sequences would be possible with i before j, but these are not Fortran-friendly.

       if (mod(count,4) == 1) then
          j1 = 2; j2 = ny-1; jstep = 1
          i1 = 2; i2 = nx-1; istep = 1
       elseif (mod(count,4) == 2) then
          j1 = ny-1; j2 = 2; jstep = -1
          i1 = nx-1; i2 = 2; istep = -1
       elseif (mod(count,4) == 3) then
          j1 = 2; j2 = ny-1; jstep = 1
          i1 = nx-1; i2 = 2; istep = -1
       elseif (mod(count,4) == 0) then
          j1 = ny-1; j2 = 2; jstep = -1
          i1 = 2; i2 = nx-1; istep = 1
       endif

       do j = j1, j2, jstep
          do i = i1, i2, istep
             if (phi_mask(i,j) == 1 .and. .not.known(i,j)) then

                ! In each cell, compute the min value of phi in the 8 neighbors
                phi_min8 = big_number
                do jj = -1,1
                   do ii = -1,1
                      ! If this is the centre point, ignore
                      if (ii == 0 .and. jj == 0) then
                         continue
                      else  ! check whether this neighbor has the minimum phi value
                         ip = i + ii
                         jp = j + jj
                         if (phi(ip,jp) < phi_min8) phi_min8 = phi(ip,jp)
                         if (mod(ii+jj,2) == 0) then  ! diagonal neighbor
                            epsilon = epsilon_diag
                         else   ! edge neighbor
                            epsilon = epsilon_edge
                         endif
                      endif
                   enddo
                enddo

                ! If phi_in(i,j) > phi_min8 + eps, set phi(i,j) = phi_in(i,j); mark cell as known.
                ! Else if phi(i,j) > phi_min8 + eps, set phi(i,j) = phi_min8 + eps; do not mark as known.
                ! Note: epsilon could be either epsilon_edge or epsilon_diag.

                if (phi_in(i,j) > phi_min8 + epsilon) then

                   !WHL - debug
                   if (verbose_depression .and. count >= 20) then
                      call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                      print*, ' '
                      print*, 'rank, i, j, ig, jg:', this_rank, i, j, iglobal, jglobal
                      print*, '   phi_in, phi:', phi_in(i,j), phi(i,j)
                      print*, '   phi_min8 =', phi_min8
                      print*, '   new phi = phi_in'
                   endif

                   phi(i,j) = phi_in(i,j)
                   known(i,j) = .true.
                   local_lowered = local_lowered + 1

                elseif (phi(i,j) > phi_min8 + epsilon) then

                   !WHL - debug
                   if (verbose_depression .and. count >= 20) then
                      call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                      print*, ' '
                      print*, 'rank, i, j, ig, jg:', this_rank, i, j, iglobal, jglobal
                      print*, '   phi_in, phi:', phi_in(i,j), phi(i,j)
                      print*, '   phi_min8 =', phi_min8
                      print*, '   new phi = phi_min8'
                   endif

                   phi(i,j) = phi_min8 + epsilon
                   local_lowered = local_lowered + 1

                endif  ! phi_in > phi_min8 + eps, phi > phi_min8 + eps

             end if   ! phi_mask = 1 and .not.known
          end do   ! i
       end do   ! j

       if (verbose_depression) then
          call point_diag(phi, 'New phi', itest, jtest, rtest, 7, 7)
       endif

       ! If one or more cells was lowered, then repeat; else exit the local loop.

       global_lowered = parallel_reduce_sum(local_lowered)

       if (global_lowered == 0) then
          finished = .true.
          if (verbose_depression .and. this_rank == rtest) then
             write(6,*) 'finished lowering'
          endif
       else
          finished = .false.
          if (verbose_depression .and. this_rank == rtest) then
             write(6,*) 'cells lowered on this iteration:', global_lowered
          endif
          call parallel_halo(phi, parallel)
       endif

       if (count > count_max) then
          call write_log('Hydrology error, exceeded max number of global iterations', GM_FATAL)
       endif

    enddo  ! finished

    if (verbose_bwat .and. this_rank == rtest) then
       write(6,*) 'Filled depressions, count =', count
    endif

  end subroutine fill_depressions

!==============================================================

  subroutine sort_heights(&
       nx,    ny,    nlocal,  &
       itest, jtest, rtest,   &
       phi,   sorted_ij)

    ! Create an array with the x and y location of each cell, sorted from from low to high values of head.
    ! Note: This subroutine sorts locally owned cells and excludes halo cells.

    ! Input/output arguments

    integer, intent(in) ::  &
         nx, ny,               & ! number of grid cells in each direction
         nlocal,               & ! number of locally owned cells
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), dimension(nx,ny), intent(in) :: &
         phi                     ! input field, to be sorted from low to high

    integer, dimension(nlocal,2), intent(inout) :: &
         sorted_ij               ! i and j indices of each cell, sorted from from low phi to high phi

    ! Local variables

    integer :: i, j, k
    integer :: ilo, ihi, jlo, jhi
    integer :: nx_local, ny_local

    real(dp), dimension(nlocal) :: vect
    integer, dimension(nlocal) :: ind

    ! Set array bounds for locally owned cells
    ilo = nhalo+1
    ihi = nx - nhalo
    jlo = nhalo+1
    jhi = ny - nhalo
    nx_local = ihi-ilo+1
    ny_local = jhi-jlo+1

    ! Fill a work vector with head values of locally owned cells
    k = 1
    do i = ilo, ihi
       do j = jlo, jhi
          vect(k) = phi(i,j)
          k = k + 1
       enddo
    enddo

    ! Sort the vector from low to high values
    ! The resulting 'ind' vector contains the k index for each cell, arranged from lowest to highest.
    ! E.g., if the lowest-ranking cell has k = 5 and the highest-ranking cell has k = 50,
    !  then ind(1) = 5 and ind(nlocal) = 50.
    ! Note: For a large problem with a small number of processors, the code can fail here
    !       because of too much recursion.

    call indexx(vect, ind)

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'Sort from low to high, nlocal =', nlocal
       print*, 'k, local i and j, ind(k), phi:'
       do k = nlocal, nlocal-10, -1
          i = floor(real(ind(k)-1)/real(ny_local)) + 1 + nhalo
          j = mod(ind(k)-1,ny_local) + 1 + nhalo
          print*, k, i, j, ind(k), phi(i,j)
       enddo
    endif

    ! Fill the sorted_ij array with the i and j values of each cell.
    ! Note: These are the i and j values we would have if there were no halo cells.
    do k = 1, nlocal
       sorted_ij(k,1) = floor(real(ind(k)-1)/real(ny_local)) + 1
       sorted_ij(k,2) = mod(ind(k)-1,ny_local) + 1
    enddo

    ! Correct the i and j values in the sorted array for halo offsets
    sorted_ij(:,:) = sorted_ij(:,:) + nhalo

  end subroutine sort_heights

!==============================================================

  subroutine get_flux_fraction(&
       nx,    ny,    nlocal,  &
       dx,    dy,             &
       itest, jtest, rtest,   &
       flux_routing_scheme,   &
       sorted_ij,             &
       head,                  &
       bwat_mask,             &
       flux_fraction)

    ! For each cell, compute the flux fraction sent to each of the 8 neighbors,
    ! based on the chosen flux routing scheme (D8, Dinf or FD8).

    ! Input/output arguments

    integer, intent(in) ::  &
         nx, ny,               & ! number of grid cells in each direction
         nlocal,               & ! number of locally owned cells
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), intent(in) ::  &
         dx, dy                  ! grid spacing in each direction (m)

    integer, intent(in) ::  &
         flux_routing_scheme     ! flux routing scheme: D8, Dinf or FD8
                                 ! D8: Flow is downhill toward the single cell with the lowest elevation.
                                 ! Dinf: Flow is downhill toward the two cells with the lowest elevations.
                                 ! FD8: Flow is downhill toward all cells with lower elevation.
                                 ! D8 scheme gives the narrowest flow, and FD8 gives the most diffuse flow.

    integer, dimension(nlocal,2), intent(in) :: &
         sorted_ij               ! i and j indices of each cell, sorted from from low phi to high phi

    real(dp), dimension(nx,ny), intent(in) :: &
         head                    ! hydraulic head (m)

    integer, dimension(nx,ny), intent(in) :: &
         bwat_mask               ! = 1 for cells in the region where basal water fluxes can be nonzero

    real(dp), dimension(-1:1,-1:1,nx,ny), intent(out) :: &
         flux_fraction           ! fraction of flux from a cell that flows downhill to each of 8 neighbors

    ! Local variables

    integer :: i, j, k, ii, jj, ip, jp, i1, i2, j1, j2, itmp, jtmp

    real(dp), dimension(-1:1,-1:1) ::  &
         dists,         &  ! distance (m) to adjacent grid cell
         slope             ! slope of head between adjacent grid cells, positive downward

    real(dp) ::  &
         slope1,        &  ! largest value of slope array
         slope2,        &  ! second largest value of slope array
         sum_slope,     &  ! sum of positive downward slopes
         slope_tmp         ! temporary slope value

    !WHL - debug
    real(dp) :: sum_frac

    ! Compute distances to adjacent grid cells for slope determination

    dists(-1,:) = (/ sqrt(dx**2 + dy**2), dy, sqrt(dx**2 + dy**2) /)
    dists(0,:) = (/ dx, 0.0d0, dx /)
    dists(1,:) = dists(-1,:)

    ! Loop through locally owned cells and compute the flux fraction sent to each neighbor cell.
    ! This fraction is stored in an array of dimension (-1:1,-1:1,nx,ny).
    ! The (0,0) element refers to the cell itself and is equal to 0 for each i and j.

    flux_fraction = 0.0d0

    do k = nlocal, 1, -1

       ! Get i and j indices of current point
       i = sorted_ij(k,1)
       j = sorted_ij(k,2)

       if (bwat_mask(i,j) == 1) then

          ! Compute the slope between this cell and each neighbor.
          ! Slopes are defined as positive for downhill neighbors, and zero otherwise.

          slope = 0.0d0

          ! Loop over adjacent points and calculate slope
          do jj = -1,1
             do ii = -1,1
                ! If this is the centre point, ignore
                if (ii == 0 .and. jj == 0) then
                   continue
                else  ! compute slope
                   ip = i + ii
                   jp = j + jj
                   if (ip >= 1 .and. ip <= nx .and. jp > 1 .and. jp <= ny) then
                      if (head(ip,jp) < head(i,j)) then
                         slope(ii,jj) = (head(i,j) - head(ip,jp)) / dists(ii,jj)
                      endif
                   endif
                endif
             enddo
          enddo

          sum_slope = sum(slope)

          if (verbose_bwat .and. this_rank == rtest .and. i == itest .and. j == jtest) then
             print*, ' '
             print*, 'slope: task, i, j =', rtest, i, j
             print*, slope(:,1)
             print*, slope(:,0)
             print*, slope(:,-1)
             print*, 'sum(slope) =', sum(slope)
          endif

          ! Distribute the downslope flux according to the flux-routing scheme:
          !  to the lowest-elevation neighbor (D8), the two lowest-elevation neighbors (Dinf), or
          !  all lower-elevation neighbors (FD8).
          ! The D8 and FD8 schemes have been tested with a simple dome problem.
          ! Dinf is less suited for the dome problem because there are many ties for 2nd greatest slope,
          !  so i2 and j2 for slope2 are not well defined.

          if (flux_routing_scheme == HO_FLUX_ROUTING_D8) then

             ! route to the adjacent cell with the lowest elevation
             slope1 = 0.0d0
             if (sum_slope > 0.d0) then
                i1 = 0; j1 = 0
                do jj = -1,1
                   do ii = -1,1
                      ip = i + ii
                      jp = j + jj
                      if (slope(ii,jj) > slope1) then
                         slope1 = slope(ii,jj)
                         i1 = ip
                         j1 = jp
                      endif
                   enddo
                enddo
             endif

             if (slope1 > 0.0d0) then
                ii = i1 - i
                jj = j1 - j
                flux_fraction(ii,jj,i,j) = 1.0d0   ! route the entire flux to one downhill cell
             else
                ! Do a fatal abort?
                print*, 'Warning: Cell with no downhill neighbors, i, j =', i, j
             endif

             if (verbose_bwat .and. this_rank == rtest .and. i == itest .and. j == jtest) then
                print*, 'i1, j1, slope1 =', i1, j1, slope1
             endif

          elseif (flux_routing_scheme == HO_FLUX_ROUTING_DINF) then

             ! route to the two adjacent cells with the lowest elevation
             i1 = 0; j1 = 0
             i2 = 0; j2 = 0
             slope1 = 0.0d0
             slope2 = 0.0d0
             do jj = -1,1
                do ii = -1,1
                   ip = i + ii
                   jp = j + jj
                   if (slope(ii,jj) > slope1) then
                      slope_tmp = slope1
                      itmp = i1
                      jtmp = j1
                      slope1 = slope(ii,jj)
                      i1 = ip
                      j1 = jp
                      slope2 = slope_tmp
                      i2 = itmp
                      j2 = jtmp
                   elseif (slope(ii,jj) > slope2) then
                      slope2 = slope(ii,jj)
                      i2 = ip
                      j2 = jp
                   endif
                enddo
             enddo

             sum_slope = slope1 + slope2    ! divide the flux between cells (i1,j1) and (i2,j2)
             if (sum_slope > 0.0d0) then
                if (slope1 > 0.0d0) then
                   ii = i1 - i
                   jj = j1 - j
                   flux_fraction(ii,jj,i,j) = slope1/sum_slope
                endif
                if (slope2 > 0.0d0) then
                   ii = i2 - i
                   jj = j2 - j
                   flux_fraction(ii,jj,i,j) = slope2/sum_slope
                endif
             else
                print*, 'Warning: Cell with no downhill neighbors, i, j =', i, j
             endif

             if (verbose_bwat .and. this_rank == rtest .and. i == itest .and. j == jtest) then
                print*, 'i1, j1, slope1:', i1, j1, slope1
                print*, 'i2, j2, slope2:', i2, j2, slope2
                print*, 'sum_slope:', sum_slope
                print*, 'slope(:, 1):', slope(:, 1)
                print*, 'slope(:, 0):', slope(:, 0)
                print*, 'slope(:,-1):', slope(:,-1)
                print*, 'flux_fraction(:, 1,i,j):', flux_fraction(:, 1,i,j)
                print*, 'flux_fraction(:, 0,i,j):', flux_fraction(:, 0,i,j)
                print*, 'flux_fraction(:,-1,i,j):', flux_fraction(:,-1,i,j)
             endif

             !WHL - bug check - make sure fractions add to 1
             sum_frac = 0.0d0
             do jj = -1,1
                do ii = -1,1
                   sum_frac = sum_frac + flux_fraction(ii,jj,i,j)
                enddo
             enddo
             if (abs(sum_frac - 1.0d0) > eps11) then
!!                print*, 'sum_frac error: r, i, j, sum:', this_rank, i, j, sum_frac
             endif

          elseif (flux_routing_scheme == HO_FLUX_ROUTING_FD8) then

             ! route to all adjacent downhill cells in proportion to grad(head)
             if (sum_slope > 0.d0) then
                do jj = -1,1
                   do ii = -1,1
                      ip = i + ii
                      jp = j + jj
                      if (slope(ii,jj) > 0.d0) then
                         flux_fraction(ii,jj,i,j) = slope(ii,jj)/sum_slope
                      endif
                   enddo
                enddo
             endif  ! sum(slope) > 0

             if (verbose_bwat .and. this_rank == rtest .and. i == itest .and. j == jtest) then
                print*, 'i1, j1, slope1:', i1, j1, slope1
                print*, 'i2, j2, slope2:', i2, j2, slope2
                print*, 'sum_slope:', sum_slope
                print*, 'slope(:, 1):', slope(:, 1)
                print*, 'slope(:, 0):', slope(:, 0)
                print*, 'slope(:,-1):', slope(:,-1)
                print*, 'flux_fraction(:, 1,i,j):', flux_fraction(:, 1,i,j)
                print*, 'flux_fraction(:, 0,i,j):', flux_fraction(:, 0,i,j)
                print*, 'flux_fraction(:,-1,i,j):', flux_fraction(:,-1,i,j)
             endif
          endif   ! flux_routing_scheme: D8, Dinf, FD8

       endif  ! bwat_mask = 1

    enddo  ! loop from high to low

  end subroutine get_flux_fraction

!==============================================================

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  ! The following two subroutines perform an index-sort of an array.
  ! They are a GPL-licenced replacement for the Numerical Recipes routine indexx.
  ! They are not derived from any NR code, but are based on a quicksort routine by
  ! Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
  ! in C, and issued under the GNU General Public License. The conversion to
  ! Fortran 90, and modification to do an index sort was done by Ian Rutt.
  !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine indexx(array, index)

    !> Performs an index sort of \texttt{array} and returns the result in
    !> \texttt{index}. The order of elements in \texttt{array} is unchanged.
    !>
    !> This is a GPL-licenced replacement for the Numerical Recipes routine indexx.
    !> It is not derived from any NR code, but are based on a quicksort routine by
    !> Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    !> in C, and issued under the GNU General Public License. The conversion to
    !> Fortran 90, and modification to do an index sort was done by Ian Rutt.

    real(dp),dimension(:) :: array !> Array to be indexed.
    integer, dimension(:) :: index !> Index of elements of \texttt{array}.
    integer :: i

    if (size(array) /= size(index)) then
       call write_log('ERROR: INDEXX size mismatch.',GM_FATAL,__FILE__,__LINE__)
    endif

    do i = 1,size(index)
       index(i) = i
    enddo

    call q_sort_index(array, index, 1, size(array))

  end subroutine indexx

!==============================================================

  recursive subroutine q_sort_index(numbers, index, left, right)

    !> This is the recursive subroutine actually used by \texttt{indexx}.
    !>
    !> This is a GPL-licenced replacement for the Numerical Recipes routine indexx.
    !> It is not derived from any NR code, but are based on a quicksort routine by
    !> Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
    !> in C, and issued under the GNU General Public License. The conversion to
    !> Fortran 90, and modification to do an index sort was done by Ian Rutt.
    !>
    !> Note: For a large problem with a small number of processors, the code can
    !        fail here with a seg fault because there is too much recursion.

    implicit none

    real(dp),dimension(:) :: numbers !> Numbers being sorted
    integer, dimension(:) :: index   !> Returned index
    integer :: left, right           !> Limit of sort region

    integer :: ll, rr
    integer :: pv_int, l_hold, r_hold, pivpos
    real(dp) :: pivot

    ll = left
    rr = right

    l_hold = ll
    r_hold = rr
    pivot = numbers(index(ll))
    pivpos = index(ll)

    do
       if (.not.(ll < rr)) exit

       do
          if  (.not.((numbers(index(rr)) >= pivot) .and. (ll < rr))) exit
          rr = rr - 1
       enddo

       if (ll /= rr) then
          index(ll) = index(rr)
          ll = ll + 1
       endif

       do
          if (.not.((numbers(index(ll)) <= pivot) .and. (ll < rr))) exit
          ll = ll + 1
       enddo

       if (ll /= rr) then
          index(rr) = index(ll)
          rr = rr - 1
       endif
    enddo

    index(ll) = pivpos
    pv_int = ll
    ll = l_hold
    rr = r_hold
    if (ll < pv_int)  call q_sort_index(numbers, index,ll, pv_int-1)
    if (rr > pv_int)  call q_sort_index(numbers, index,pv_int+1, rr)

  end subroutine q_sort_index

!==============================================================

end module glissade_basal_water

!==============================================================
