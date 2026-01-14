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

   use glimmer_global, only: dp, i8
   use glimmer_paramets, only: iulog, eps11, eps08
   use glimmer_physcon, only: rhoi, rhow, lhci, grav, scyr
   use glimmer_log
   use glimmer_utils, only: point_diag
   use glide_types
   use cism_parallel, only: main_task, this_rank, nhalo, parallel_type, &
        parallel_halo, parallel_global_sum

   !WHL - debug
   use glimmer_utils, only: double_to_binary

   implicit none

   private
   public :: glissade_basal_water_init, glissade_calcbwat, glissade_bwat_flux_routing

!!   logical, parameter :: verbose_bwat = .false.
   logical, parameter :: verbose_bwat = .true.

   character(len=64) :: binary_str

   ! two versions of this subroutine; the second supports reproducible sums
   interface route_flux_to_margin_or_halo
      module procedure route_flux_to_margin_or_halo_real8
      module procedure route_flux_to_margin_or_halo_integer8
   end interface

 contains

!==============================================================

  subroutine glissade_basal_water_init(model)

    type(glide_global_type) :: model

    select case (model%options%which_ho_bwat)

    ! HO_BWAT_NONE:         basal water depth = 0
    ! HO_BWAT_CONSTANT:     basal water depth = prescribed constant
    ! HO_BWAT_LOCAL_TILL:   local basal till model with prescribed drainage rate
    ! HO_BWAT_FLUX_ROUTING: steady-state water routing with flux calculation

    case(HO_BWAT_CONSTANT)

       ! Set a constant water thickness where ice is present
       where (model%geometry%thck > model%numerics%thklim)
          model%basal_hydro%bwat(:,:) = model%basal_hydro%const_bwat
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
       delta_Tb,                     &
       btemp_flow_scale,             &
       btemp_freeze_scale,           &
       bwatflx,       bwat_diag,     &
       bhydroflx,                    &
       head,          grad_head,     &
       reprosum_in)

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
         bmlt_hydro,              & ! meltwater input rate (m/s); can include an englacial or surface source
         delta_Tb                   ! difference T_pmp - T_bed (degC)

    real(dp), intent(in) ::  &
         thklim                     ! minimum ice thickness for basal melt and hydropotential calculations (m)
                                    ! Note: This is typically model%geometry%thklim_temp

    ! Note: These scales ensure a smooth transition in behavior between frozen and thawed beds.
    !       Both scales are computed in a similar way, but they apply to different parts of the algorithm.
    ! TODO: Decide whether to keep both scales. Only the flow scale works for reprosums, so we might want
    !       to remove the freeze scale.
    real(dp), intent(in) ::  &
         btemp_flow_scale,        & ! temperature scale for routing water flow around cells with a frozen bed (deg C)
         btemp_freeze_scale         ! temperature scale for refreezing water beneath cells with a frozen bed (degC)

    integer, dimension(nx,ny), intent(in) ::  &
         bwat_mask,               & ! mask to identify cells through which basal water is routed;
                                    ! = 0 for floating and ocean cells, cells at global domain edge,
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
         bwatflx,                 & ! basal water flux through each grid cell (m/s)
         bwat_diag,               & ! diagnosed basal water depth (m)
         bhydroflx,               & ! basal heat flux from refreezing meltwater (W/m2); enters thermal solve later
         head,                    & ! hydraulic head (m)
         grad_head                  ! gradient of hydraulic head (m/m), averaged to cell centers

    ! Note: The reprosum option requires (1) D8 routing (each cell routes its flux to one downstream neighbor only)
    !       and (2) no temperature-weighted refreezing.
    logical, intent(in), optional :: &
         reprosum_in                ! if true, then do a computation independent of the number of tasks

    ! Local variables

    integer :: i, j, p

    real(dp), dimension(nx,ny) :: &
         bwatflx_refreeze,        & ! water flux held for refreezing (m^3/s)
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

    logical :: reprosum                      ! local version of reprosum_in

    integer :: nx_test, ny_test
    real(dp), dimension(:,:), allocatable :: phi_test
    integer,  dimension(:,:), allocatable :: mask_test

    if (verbose_bwat .and. this_rank == rtest) then
       write(iulog,*) 'In glissade_bwat_flux_routing: rtest, itest, jtest =', rtest, itest, jtest
    endif

    if (present(reprosum_in)) then
       reprosum = reprosum_in
    else
       reprosum = .false.
    endif

    ! Uncomment if the following fields are not already up to date in halo cells
!    call parallel_halo(thk,  parallel)
!    call parallel_halo(topg, parallel)
    call parallel_halo(bmlt_hydro, parallel)

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

    call route_basal_water(&
         nx,      ny,            &
         dx,      dy,            &
         parallel,               &
         itest, jtest, rtest,    &
         flux_routing_scheme,    &
         head,                   &
         bmlt_hydro,             &
         delta_Tb,               &
         btemp_flow_scale,       &
         btemp_freeze_scale,     &
         bwat_mask,              &
         bwatflx,                &
         bwatflx_refreeze,       &
         lakes,                  &
         reprosum)

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

    ! Convert bwatflx from m^3/s to m/s for output
    bwatflx(:,:) = bwatflx(:,:) / (dx*dy)

    ! Given bwatflx_refreeze in m^3/s, compute bhydroflx in W/m2.
    ! This is the heat flux needed to refreeze the meltwater held in each cell.
    ! This heat flux is supplied at the bed during the next thermal solve.
    ! If there is more than enough heat to thaw the bed, some meltwater will be returned later
    !  instead of refrozen.

    bhydroflx(:,:) = bwatflx_refreeze(:,:) * rhoi * lhci / (dx*dy)

    if (verbose_bwat) then
       call point_diag(bwatflx*scyr, 'Final bwatflx (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(bwatflx_refreeze*scyr/(dx*dy), 'bwatflx_refreeze (m/yr)', itest, jtest, rtest, 7, 7)
       call point_diag(bhydroflx, 'bhydroflx (W/m2)', itest, jtest, rtest, 7, 7, '(f10.6)')
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
         delta_Tb,               &
         btemp_flow_scale,       &
         btemp_freeze_scale,     &
         bwat_mask,              &
         bwatflx,                &
         bwatflx_refreeze,       &
         lakes,                  &
         reprosum_in)

    ! Route water from the basal melt field to its destination, recording the water flux along the way.
    ! Water flow direction is determined according to the gradient of the hydraulic head.
    ! For the algorithm to work correctly, surface depressions must be filled,
    !  so that all cells have an outlet to the ice sheet margin.
    ! This results in the lakes field, which is the difference between the filled head and the original head.
    !  The method used is by Quinn et. al. (1991).
    !
    ! Originally based on code by Jesse Johnson and Ian Rutt in the Glimmer model

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
         bmlt_hydro,           & ! meltwater input rate (m/s)
         delta_Tb                ! difference T_pmp - T_bed (degC)

    real(dp), intent(in) :: &
         btemp_flow_scale,     & ! temperature scale for routing water flow around cells with a frozen bed (deg C)
         btemp_freeze_scale      ! temperature scale for refreezing water beneath cells with a frozen bed (degC)
                                 ! If scale = 0, assume no temperature dependence

    real(dp), dimension(nx,ny), intent(inout)  ::  &
         head                    ! hydraulic head (m)
                                 ! intent inout because it can be filled and adjusted below

    integer, dimension(nx,ny), intent(in)  ::  &
         bwat_mask               ! mask to identify cells through which basal water is routed;
                                 ! excludes floating and ocean cells

    real(dp), dimension(nx,ny), intent(out) ::  &
         bwatflx,             &  ! water flux through a grid cell (m^3/s)
         bwatflx_refreeze,    &  ! water flux held for refreezing (m^3/s)
         lakes                   ! lakes field, difference between filled and original head

    logical, intent(in), optional :: &
         reprosum_in             ! if true, then do a computation independent of the number of tasks

    ! Local variables

    integer :: nlocal            ! number of locally owned cells
    integer :: count, count_max  ! iteration counters
    integer :: i, j, k, iglobal, jglobal
    integer :: ii, jj, imax, jmax

    logical :: finished    ! true when an iterative loop has finished

    integer,  dimension(:,:), allocatable ::  &
         sorted_ij         ! i and j indices of all cells, sorted from low to high values of head

    real(dp), dimension(-1:1,-1:1,nx,ny) ::  &
         flux_fraction     ! fraction of flux from each cell that flows downhill to each of 8 neighbors

    real(dp), dimension(nx,ny) ::  &
         head_filled,           & ! head after depressions are filled (m)
         btemp_weight_flow,     & ! temp-dependent weighting factor, forcing flow around cells with frozen beds
         btemp_weight_freeze,   & ! temp-dependent weighting factor, favoring refreezing in cells with frozen beds
         bwatflx_accum,         & ! water flux through the cell (m^3/s) accumulated over multiple iterations
         bwatflx_refreeze_accum   ! water flux (m^3/s) refreezing in place, accumulated over multiple iterations

    integer, dimension(nx,ny) ::  &
         local_mask,            & ! = 1 for cells owned by the local processor, else = 0
         halo_mask,             & ! = 1 for the layer of halo cells adjacent to locally owned cells, else = 0
         margin_mask              ! = 1 for cells at the grounded ice margin, as defined by bwat_mask, else = 0

    real(dp) :: &
         total_flux_in,         & ! total input flux (m^3/s), computed as sum of bmlt_hydro*dx*dy
         total_flux_margin,     & ! total output flux (m^3/s) at the ice margin
         total_flux_refreeze,   & ! total flux (m^3/s) to refreeze internally
         total_flux_out,        & ! sum of total_bwatflx_margin and total_bwatflx_refreeze
         err,                   & ! water conservation error
         global_flux_sum          ! flux sum over all cells in global domain

    ! The following i8 variables are for computing reproducible sums
    integer(i8), dimension(nx,ny) ::  &
         bwatflx_int,               & ! water flux through a grid cell (m^3/s)
         btemp_weight_freeze_int,   & ! temp-dependent weighting factor, favoring refreezing in cells with frozen beds
         bwatflx_accum_int,         & ! water flux through the cell (m^3/s) accumulated over multiple iterations
         bwatflx_refreeze_accum_int   ! water flux (m^3/s) refreezing in place, accumulated over multiple iterations

    integer(i8), dimension(-1:1,-1:1,nx,ny) ::  &
         flux_fraction_int            ! fraction of flux from each cell that flows downhill to each of 8 neighbors

    real(dp), parameter :: &
         factor_bwatflx = 1.d16       ! factor for converting between bwatflx and bwatflx_int;
                                      ! large value desired for water mass conservation

    logical :: reprosum               ! local version of reprosum_in

    character(len=100) :: message

    !WHL - debug
    character(len=64) :: binary_str

    if (present(reprosum_in)) then
       reprosum = reprosum_in
    else
       reprosum = .false.
    endif

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
       write(iulog,*) ' '
       write(iulog,*) 'sorted, from the top:'
       do k = nlocal, nlocal-10, -1
          i = sorted_ij(k,1)
          j = sorted_ij(k,2)
          write(iulog,*) k, i, j, head(i,j)
       enddo
    endif

    ! Make sure bmlt_hydro >= 0 everywhere
    do j = 1, ny
       do i = 1, nx
          if (bmlt_hydro(i,j) < 0.0d0) then
             call parallel_globalindex(i, j, iglobal, jglobal, parallel)
!             write(iulog,*) 'Hydrology error: bmlt_hydro < 0, iglobal, jglobal, bmlt_hydro =', &
!                     iglobal, jglobal, bmlt_hydro(i,j)
             write(message,*) 'Hydrology error: bmlt_hydro < 0, iglobal, jglobal, bmlt_hydro =', &
                  iglobal, jglobal, bmlt_hydro(i,j)
             call write_log(message, GM_FATAL)
          endif
       enddo
    enddo

    ! Compute temperature-dependent weighting factors for flux routing.
    ! There are two scales with related but distinct functions:
    ! (1) In subroutine get_flux_fraction, btemp_flow_scale is used to weigh potential downstream paths.
    !     A low value of btemp_weight_flow means that water is less likely to pass through.
    ! (2) When water enters a frozen cell (delta_Tb > 0), btemp_freeze_scale determines
    !     how much of the flux refreezes in place rather than passing through.
    !     A small value of btemp_weight_freeze means that more water refreezes, and less passes through.
    ! Note: For reproducible sums, refreezing is not supported; must have btemp_weight_freeze = 1 everywhere.
    ! TODO: Possibly remove btemp_freeze_scale and just keep btemp_flow_scale.

    btemp_weight_flow = 1.0d0
    if (btemp_flow_scale > 0.0d0) then
       if (.not. reprosum) then
          where (bwat_mask == 1)
             where (delta_Tb > 0.0d0)
                btemp_weight_flow = exp(-delta_Tb/btemp_flow_scale)
             endwhere
          endwhere
       endif
    endif

    btemp_weight_freeze = 1.0d0
    if (btemp_freeze_scale > 0.0d0) then
       if (.not. reprosum) then
          where (bwat_mask == 1)
             where (delta_Tb > 0.0d0)
                btemp_weight_freeze = exp(-delta_Tb/btemp_freeze_scale)
             endwhere
          endwhere
       endif
    endif

    if (verbose_bwat) then
       call point_diag(delta_Tb, 'Tpmp - Tb', itest, jtest, rtest, 7, 7)
       call point_diag(btemp_weight_flow, 'btemp_weight_flow', itest, jtest, rtest, 7, 7)
       call point_diag(btemp_weight_freeze, 'btemp_weight_freeze', itest, jtest, rtest, 7, 7)
    endif

    ! Compute the fraction of the incoming flux sent to each downstream neighbor.

    call get_flux_fraction(&
         nx,    ny,    nlocal,  &
         dx,    dy,             &
         itest, jtest, rtest,   &
         flux_routing_scheme,   &
         sorted_ij,             &
         head,                  &
         btemp_weight_flow,     &
         bwat_mask,             &
         flux_fraction)

    ! Initialize bwatflx in locally owned cells.
    ! Set to the local melt rate, multiplied by area (so the units are m^3/s).

    bwatflx = 0.0d0
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          bwatflx(i,j) = bmlt_hydro(i,j) * dx * dy
       enddo
    enddo

    ! Compute total input of meltwater (m^3/s)
    total_flux_in = parallel_global_sum(bwatflx, parallel)
    if (verbose_bwat .and. this_rank == rtest) then
       write(iulog,*) 'Total input basal melt flux (m^3/s):', total_flux_in
!!       call double_to_binary(total_flux_in, binary_str)
!!       write(iulog,*) '   Binary string', binary_str
    endif

    ! Route the water downstream, keeping track of the steady-state flux through each cell.
    ! The loop goes from highest to lowest values of head on the local processor.
    ! At the end of the loop, all the incoming flux has either been
    ! (1) routed to the ice sheet margin,
    ! (2) set aside for later refreezing, or
    ! (3) routed to a halo cell, from which it will continue downstream on the next iteration.
    ! When all the water has been routed to the margin or set aside for refreezing, we are done.

    ! Note: It is hard to predict how many iterations will be sufficient.
    !       With Dinf or FD8 we can have flow back and forth across processor boundaries,
    !        requiring many iterations to reach the margin.
    !       For Greenland 4 km, Dinf requires ~20 iterations on 4 cores, and FD8 can require > 40.
    !       For Antarctica 8 km, FD8 can require > 50.

    ! Initialize the cumulative fluxes
    bwatflx_accum = 0.0d0
    bwatflx_refreeze_accum = 0.0d0

    if (reprosum) then

       ! Convert bwatflx to a scaled i8 array
       bwatflx_int(:,:) = nint(bwatflx(:,:)*factor_bwatflx, i8)

       ! Convert flux_fraction to i8
       ! Note: This will work only for the D8 scheme, where all the flux goes downstream
       !       to a single cell.
       flux_fraction_int(:,:,:,:) = nint(flux_fraction(:,:,:,:), i8)

       ! Convert btemp_weight_freeze to i8
       btemp_weight_freeze_int(:,:) = 1
       ! Note: Can round up to 1 and down to 0 by uncommenting the following line.
       !       However, a mix of 1's and 0's leads to oscillations in basal temperature,
       !        so it is safer to turn off refreezing by setting btemp_weight_freeze = 1 everywhere.
!       btemp_weight_freeze_int(:,:) = nint(btemp_weight_freeze(:,:), i8)

       ! Initialize other arrays
       bwatflx_accum_int = 0
       bwatflx_refreeze_accum_int = 0

    endif   ! reprosum

    count = 0
    count_max = 100
    finished = .false.

    do while (.not.finished)

       count = count + 1
       if (verbose_bwat .and. this_rank == rtest) write(iulog,*) 'flux routing, count =', count
       if (count > count_max) then
          call write_log('Hydrology error: too many iterations in route_basal_water', GM_FATAL)
       endif

       if (reprosum) then

          ! route downstream
          ! Note: The fluxes are scaled by factor_bwatflx

          call route_flux_to_margin_or_halo(&
            nx, ny, nlocal,             &
            itest, jtest, rtest, count, &
            parallel,                   &
            sorted_ij,                  &
            local_mask,                 &
            halo_mask,                  &
            bwat_mask,                  &
            flux_fraction_int,          &
            btemp_weight_freeze_int,    &
            bwatflx_int,                &
            bwatflx_accum_int,          &
            bwatflx_refreeze_accum_int, &
            finished)

          if (verbose_bwat .and. this_rank == rtest .and. count <= 2) then
             i = itest; j = jtest
             write(iulog,*) 'count, rank i, j, bwatflx_accum (m/yr), bwatflx_refreeze_accum:', &
                  count, rtest, i, j, real(bwatflx_accum_int(i,j),dp)/factor_bwatflx, &
                  real(bwatflx_refreeze_accum_int(i,j),dp)/factor_bwatflx
          endif

       else   ! non-reproducible sums

          call route_flux_to_margin_or_halo(&
            nx, ny, nlocal,      &
            itest, jtest, rtest, count, &
            parallel,            &
            sorted_ij,           &
            local_mask,          &
            halo_mask,           &
            bwat_mask,           &
            flux_fraction,       &
            btemp_weight_freeze, &
            bwatflx,             &
            bwatflx_accum,       &
            bwatflx_refreeze_accum, &
            finished)

          if (verbose_bwat .and. this_rank == rtest .and. count <= 2) then
             i = itest; j = jtest
             write(iulog,*) 'count, rank i, j, bwatflx_accum(m/yr), bwatflx_refreeze_accum:', &
                  count, rtest, i, j, bwatflx_accum(i,j) * scyr/(dx*dy), &
                  bwatflx_refreeze_accum(i,j) * scyr/(dx*dy)
          endif

       endif   ! reprosum

    enddo   ! finished

    if (reprosum) then
       ! Convert fluxes back to real(dp)
       bwatflx_accum = real(bwatflx_accum_int, dp) / factor_bwatflx
       bwatflx_refreeze_accum = real(bwatflx_refreeze_accum_int, dp) / factor_bwatflx
    endif

    ! Copy the accumulated values to the output arrays bwatflx and bwatflx_refreeze
    bwatflx = bwatflx_accum
    bwatflx_refreeze = bwatflx_refreeze_accum
    if (verbose_bwat .and. this_rank == rtest) write(iulog,*) 'Done routing fluxes'

    ! Identify cells just beyond the ice sheet margin, which can receive from upstream but not send downstream
    where (bwat_mask == 0 .and. bwatflx > 0.0d0)
       margin_mask = 1
    elsewhere
       margin_mask = 0
    endwhere

    ! Compute total output of meltwater (m^3/s) and check that input = output, within roundoff.
    total_flux_margin   = parallel_global_sum(bwatflx*margin_mask, parallel)
    total_flux_refreeze = parallel_global_sum(bwatflx_refreeze*(1.0d0 - margin_mask), parallel)
    total_flux_out = total_flux_margin + total_flux_refreeze

    if (verbose_bwat .and. this_rank == rtest) then
       write(iulog,*) 'Total bwatflx at margin (m^3/s):', total_flux_margin
       write(iulog,*) 'Total bwatflx_refreeze (m^3/s)=', total_flux_refreeze
       write(iulog,*) 'Total bwatflx (m^3/s)=', total_flux_out
       write(iulog,*) 'Difference between output and input =', total_flux_out - total_flux_in
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
         nx,       ny,        &
         dx,       dy,        &
         itest, jtest, rtest, &
         head,                &
         dhead_dx, dhead_dy,  &
         ice_mask,            &
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

    use cism_parallel, only: parallel_reduce_sum, parallel_globalindex

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

    ! Local variables

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

    logical, parameter :: verbose_depression = .false.

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
          write(iulog,*) ' '
          write(iulog,*) 'fill_depressions, count =', count
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
                      write(iulog,*) ' '
                      write(iulog,*) 'rank, i, j, ig, jg:', this_rank, i, j, iglobal, jglobal
                      write(iulog,*) '   phi_in, phi:', phi_in(i,j), phi(i,j)
                      write(iulog,*) '   phi_min8 =', phi_min8
                      write(iulog,*) '   new phi = phi_in'
                   endif

                   phi(i,j) = phi_in(i,j)
                   known(i,j) = .true.
                   local_lowered = local_lowered + 1

                elseif (phi(i,j) > phi_min8 + epsilon) then

                   !WHL - debug
                   if (verbose_depression .and. count >= 20) then
                      call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                      write(iulog,*) ' '
                      write(iulog,*) 'rank, i, j, ig, jg:', this_rank, i, j, iglobal, jglobal
                      write(iulog,*) '   phi_in, phi:', phi_in(i,j), phi(i,j)
                      write(iulog,*) '   phi_min8 =', phi_min8
                      write(iulog,*) '   new phi = phi_min8'
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
             write(iulog,*) 'finished lowering'
          endif
       else
          finished = .false.
          if (verbose_depression .and. this_rank == rtest) then
             write(iulog,*) 'cells lowered on this iteration:', global_lowered
          endif
          call parallel_halo(phi, parallel)
       endif

       if (count > count_max) then
          call write_log('Hydrology error, exceeded max number of global iterations', GM_FATAL)
       endif

    enddo  ! finished

    if (verbose_bwat .and. this_rank == rtest) then
       write(iulog,*) 'Filled depressions, count =', count
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
         itest, jtest, rtest     ! coordinates of diagnostic point  !! not currently used

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
       btemp_weight_flow,     &
       bwat_mask,             &
       flux_fraction)

    ! For each cell, compute the flux fraction sent to each of the 8 neighbors,
    !  based on the chosen flux routing scheme (D8, Dinf or FD8).
    ! The flux fraction is proportional to grad(head).

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
         head,                 & ! hydraulic head (m)
         btemp_weight_flow       ! temperature-dependent weighting factor, forcing flow around cells with frozen beds

    integer, dimension(nx,ny), intent(in) :: &
         bwat_mask               ! = 1 for cells in the region where basal water fluxes can be nonzero

    real(dp), dimension(-1:1,-1:1,nx,ny), intent(out) :: &
         flux_fraction           ! fraction of flux from a cell that flows downhill to each of 8 neighbors

    ! Local variables

    integer :: i, j, k, ii, jj, ip, jp, i1, i2, j1, j2, itmp, jtmp

    real(dp), dimension(-1:1,-1:1) ::  &
         dists,              & ! distance (m) to adjacent grid cell
         slope                 ! slope of head between adjacent grid cells, positive downward

    real(dp) ::  &
         slope1,             & ! largest value of slope array
         slope2,             & ! second largest value of slope array
         sum_slope,          & ! sum of positive downward slopes
         slope_tmp             ! temporary slope value

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
                else  ! compute the hydropotential slope, weighted based on downstream bed temperature
                   ip = i + ii
                   jp = j + jj
                   if (ip >= 1 .and. ip <= nx .and. jp > 1 .and. jp <= ny) then
                      if (head(ip,jp) < head(i,j)) then
                         slope(ii,jj) = btemp_weight_flow(ip,jp) * (head(i,j) - head(ip,jp)) / dists(ii,jj)
                      endif
                   endif
                endif
             enddo
          enddo

          sum_slope = sum(slope)

          if (verbose_bwat .and. this_rank == rtest .and. i == itest .and. j == jtest) then
             write(iulog,*) ' '
             write(iulog,*) 'slope: task, i, j =', rtest, i, j
             write(iulog,*) slope(:,1)
             write(iulog,*) slope(:,0)
             write(iulog,*) slope(:,-1)
             write(iulog,*) 'sum(slope) =', sum(slope)
          endif

          ! Distribute the downslope flux according to the flux-routing scheme:
          !  to the lowest-elevation neighbor (D8), the two lowest-elevation neighbors (Dinf), or
          !  all lower-elevation neighbors (FD8).
          ! The D8 and FD8 schemes have been tested with a simple dome problem.
          ! Dinf is less suited for the dome problem because there are many ties for 2nd greatest slope,
          !  so i2 and j2 for slope2 are not well defined.

          if (flux_routing_scheme == HO_FLUX_ROUTING_D8) then

             ! route to the lowest-lying adjacent cell
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
                write(iulog,*) 'Warning: Cell with no downhill neighbors, i, j =', i, j
             endif

             if (verbose_bwat .and. this_rank == rtest .and. i == itest .and. j == jtest) then
                write(iulog,*) 'i1, j1, slope1 =', i1, j1, slope1
             endif

          elseif (flux_routing_scheme == HO_FLUX_ROUTING_DINF) then

             ! route to the two lowest-lying adjacent cells
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
                write(iulog,*) 'Warning: Cell with no downhill neighbors, i, j =', i, j
             endif

             if (verbose_bwat .and. this_rank == rtest .and. i == itest .and. j == jtest) then
                write(iulog,*) 'i1, j1, slope1:', i1, j1, slope1
                write(iulog,*) 'i2, j2, slope2:', i2, j2, slope2
                write(iulog,*) 'sum_slope:', sum_slope
                write(iulog,*) 'slope(:, 1):', slope(:, 1)
                write(iulog,*) 'slope(:, 0):', slope(:, 0)
                write(iulog,*) 'slope(:,-1):', slope(:,-1)
                write(iulog,*) 'flux_fraction(:, 1,i,j):', flux_fraction(:, 1,i,j)
                write(iulog,*) 'flux_fraction(:, 0,i,j):', flux_fraction(:, 0,i,j)
                write(iulog,*) 'flux_fraction(:,-1,i,j):', flux_fraction(:,-1,i,j)
             endif

             !WHL - bug check - make sure fractions add to 1
             sum_frac = 0.0d0
             do jj = -1,1
                do ii = -1,1
                   sum_frac = sum_frac + flux_fraction(ii,jj,i,j)
                enddo
             enddo
             if (abs(sum_frac - 1.0d0) > eps11) then
!!                write(iulog,*) 'sum_frac error: r, i, j, sum:', this_rank, i, j, sum_frac
             endif

          elseif (flux_routing_scheme == HO_FLUX_ROUTING_FD8) then

             ! route to all adjacent downhill cells in proportion to the slope
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
                write(iulog,*) 'i1, j1, slope1:', i1, j1, slope1
                write(iulog,*) 'i2, j2, slope2:', i2, j2, slope2
                write(iulog,*) 'sum_slope:', sum_slope
                write(iulog,*) 'slope(:, 1):', slope(:, 1)
                write(iulog,*) 'slope(:, 0):', slope(:, 0)
                write(iulog,*) 'slope(:,-1):', slope(:,-1)
                write(iulog,*) 'flux_fraction(:, 1,i,j):', flux_fraction(:, 1,i,j)
                write(iulog,*) 'flux_fraction(:, 0,i,j):', flux_fraction(:, 0,i,j)
                write(iulog,*) 'flux_fraction(:,-1,i,j):', flux_fraction(:,-1,i,j)
             endif

          endif   ! flux_routing_scheme: D8, Dinf, FD8

       endif  ! bwat_mask = 1

    enddo  ! loop from high to low

  end subroutine get_flux_fraction

!==============================================================

  subroutine route_flux_to_margin_or_halo_real8(&
       nx, ny, nlocal,      &
       itest, jtest, rtest, count, &
       parallel,            &
       sorted_ij,           &
       local_mask,          &
       halo_mask,           &
       bwat_mask,           &
       flux_fraction,       &
       btemp_weight_freeze, &
       bwatflx,             &
       bwatflx_accum,       &
       bwatflx_refreeze_accum, &
       finished)

    ! Given the input bwatflx, route the water downstream, keeping track of fluxes along the way.
    ! The loop goes from highest to lowest values of 'head' on the local processor.
    ! At the end of the loop, all the incoming flux has either been
    ! (1) routed to the ice sheet margin;
    ! (2) set aside for later refreezing; or
    ! (3) routed to a halo cell, from which it continues downstream the next time the subroutine is called.
    ! The subroutine is called iteratively until all no water remains in halo cells.

    implicit none

    ! Input/output variables

    integer, intent(in) :: &
         nx, ny,               & ! number of cells in each direction
         nlocal,               & ! number of locally owned grid cells on the processor
         itest, jtest, rtest,  & ! coordinates of diagnostic point
         count                   ! iteration count (diagnostic only)

    type(parallel_type), intent(in) ::  &
         parallel                ! info for parallel communication

    integer, dimension(nlocal,2), intent(in)  :: &
         sorted_ij               ! i and j indices of each local cell, sorted low to high

    integer, dimension(nx,ny), intent(in) :: &
         local_mask,           & ! = 1 for cells owned by the local processor, else = 0
         halo_mask,            & ! = 1 for the layer of halo cells adjacent to locally owned cells, else = 0
         bwat_mask               ! = 1 for cells through which basal water is routed; excludes floating and ocean cells

    real(dp), dimension(-1:1,-1:1,nx,ny), intent(in) ::  &
         flux_fraction           ! fraction of flux from a cell that flows downhill to each of 8 neighbors
                                 ! last two indices identify the source cell;
                                 ! 1st two indices give relative location of receiving cell

    real(dp), dimension(nx,ny), intent(in) :: &
         btemp_weight_freeze     ! temperature-dependent weighting factor, favoring refreezing at frozen beds

    real(dp), dimension(nx,ny), intent(inout) :: &
         bwatflx,              & ! on input: water flux (m^3/s) to be routed to the margin or halo
                                 ! on output: flux routed to halo, to be routed further next time
         bwatflx_accum,        & ! cumulative bwatflx (m/3/s) over multiple iterations
         bwatflx_refreeze_accum  ! cumulative bwatflx_refreeze (m^s/s) over multiple iterations

    logical, intent(inout) :: &
         finished                ! initially F; set to T when all water has been routed as far as it can go

    ! Local variables

    integer :: i, j, k
    integer :: ii, jj, ip, jp

    real(dp), dimension(-1:1,-1:1,nx,ny)::  &
         bwatflx_halo            ! flux routed to halo cells
                                 ! last two indices identify the source cell;
                                 ! 1st two indices give relative location of receiving cell

    real(dp), dimension(nx,ny) :: &
         bwatflx_refreeze,     & ! flux (m^3/s) saved for later refreezing; not routed further downstream
         sum_bwatflx_halo        ! bwatflx_halo summed over the first 2 indices

    real(dp) :: &
         flx_thru,             & ! flux (m^3/s) that continues downstream
         global_halo_sum         ! global sum of water flux in halo cells

    ! Initialize fluxes
    bwatflx_refreeze = 0.0d0
    bwatflx_halo = 0.0d0

    ! loop from high to low values on the local processor
    do k = nlocal, 1, -1

       ! Get i and j indices of current cell
       i = sorted_ij(k,1)
       j = sorted_ij(k,2)

       if (bwat_mask(i,j) == 1 .and. bwatflx(i,j) > 0.0d0) then

          ! Distribute the flux to downstream neighbors.
          ! Based on the temperature-dependent weighting factor btemp_weight_freeze, all or part of the flux
          !  is refrozen in place instead of being routed downstream.
          flx_thru = bwatflx(i,j) * btemp_weight_freeze(i,j)
          bwatflx_refreeze(i,j) = bwatflx(i,j) * (1.0d0 - btemp_weight_freeze(i,j))
          do jj = -1,1
             do ii = -1,1
                ip = i + ii
                jp = j + jj
                if (flux_fraction(ii,jj,i,j) > 0.0d0) then
                   if (halo_mask(ip,jp) == 1) then
                      bwatflx_halo(ii,jj,i,j) = flx_thru*flux_fraction(ii,jj,i,j)
                      if (verbose_bwat .and. this_rank==rtest .and. i==itest .and. j==jtest .and. count <= 2) then
                         write(iulog,*) 'Flux to halo, i, j, ii, jj, flux:', &
                              i, j, ii, jj, flx_thru*flux_fraction(ii,jj,i,j)
                      endif
                   elseif (local_mask(ip,jp) == 1) then
                      bwatflx(ip,jp) = bwatflx(ip,jp) + flx_thru*flux_fraction(ii,jj,i,j)
                      if (verbose_bwat .and. this_rank==rtest .and. i==itest .and. j==jtest .and. count <= 2) then
                         write(iulog,*) 'Flux to neighbor, i, j, ii, jj, flux:', &
                              i, j, ii, jj, flx_thru*flux_fraction(ii,jj,i,j)
                      endif
                   endif
                endif   ! flux_fraction > 0
             enddo  ! ii
          enddo  ! jj
       endif  ! bwat_mask = 1, bwatflx > 0
    enddo  ! loop from high to low

    ! Accumulate the fluxes in the output arrays
    bwatflx_accum = bwatflx_accum + bwatflx
    bwatflx_refreeze_accum = bwatflx_refreeze_accum + bwatflx_refreeze

    ! Compute the total bwatflx in halo cells
    do j = 1, ny
       do i = 1, nx
          sum_bwatflx_halo(i,j) = sum(bwatflx_halo(:,:,i,j))
       enddo
    enddo
    global_halo_sum = parallel_global_sum(sum_bwatflx_halo, parallel)

    ! If bwatflx_halo = 0 everywhere, then we are done.
    ! Where bwatflx_halo is nonzero, communicate it to the neighboring task.
    ! It will be routed further downstream the next time this subroutine is called.

    if (global_halo_sum > 0.0d0) then

       if (verbose_bwat .and. count <= 2) then
          if (this_rank == rtest) write(iulog,*) 'Before halo update, global_halo_sum:', global_halo_sum
          call point_diag(sum_bwatflx_halo, 'sum_bwatflx_halo', itest, jtest, rtest, 7, 7)
       endif

       ! Reset bwatflx to zero for the halo transfer
       bwatflx = 0.0d0

       ! Communicate bmltflx_halo to the halo cells of neighboring processors
       call parallel_halo(bwatflx_halo(:,:,:,:), parallel)

       ! bmltflx_halo is now available in the halo cells of the local processor.
       ! Route downslope to the adjacent locally owned cells.
       ! These fluxes will be routed further downstream during the next iteration.

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
                            if (verbose_bwat .and. ip==itest .and. jp==jtest .and. this_rank==rtest .and. count <= 2) then
                               write(iulog,*) 'Nonzero bwatflx from halo, rank, i, j:', &
                                    this_rank, ip, jp, bwatflx_halo(ii,jj,i,j)
                            endif
                         endif
                      endif   !  bwatflx_halo > 0 to a local cell
                   enddo   ! ii
                enddo   ! jj
             endif   ! bwatflx_halo > 0 from this halo cell
          enddo   ! i
       enddo   ! j

    else

       finished = .true.   ! no water in halo cells to route further

    endif  ! global_halo_sum > 0

  end subroutine route_flux_to_margin_or_halo_real8

!==============================================================

  subroutine route_flux_to_margin_or_halo_integer8(&
       nx, ny, nlocal,             &
       itest, jtest, rtest, count, &
       parallel,                   &
       sorted_ij,                  &
       local_mask,                 &
       halo_mask,                  &
       bwat_mask,                  &
       flux_fraction,              &
       btemp_weight_freeze,        &
       bwatflx,                    &
       bwatflx_accum,              &
       bwatflx_refreeze_accum,     &
       finished)

    ! Given the input bwatflx, route the water downstream, keeping track of fluxes along the way.
    ! The loop goes from highest to lowest values of 'head' on the local processor.
    ! At the end of the loop, all the incoming flux has either been
    ! (1) routed to the ice sheet margin;
    ! (2) set aside for later refreezing; or
    ! (3) routed to a halo cell, from which it continues downstream the next time the subroutine is called.
    ! The subroutine is called iteratively until all no water remains in halo cells.

    implicit none

    ! Input/output variables

    integer, intent(in) :: &
         nx, ny,               & ! number of cells in each direction
         nlocal,               & ! number of locally owned grid cells on the processor
         itest, jtest, rtest,  & ! coordinates of diagnostic point
         count                   ! iteration count (diagnostic only)

    type(parallel_type), intent(in) ::  &
         parallel                ! info for parallel communication

    integer, dimension(nlocal,2), intent(in)  :: &
         sorted_ij               ! i and j indices of each local cell, sorted low to high

    integer, dimension(nx,ny), intent(in) :: &
         local_mask,           & ! = 1 for cells owned by the local processor, else = 0
         halo_mask,            & ! = 1 for the layer of halo cells adjacent to locally owned cells, else = 0
         bwat_mask               ! = 1 for cells through which basal water is routed; excludes floating and ocean cells

    ! Note: Both flux_fraction and btemp_weight_freeze are constrained to be 0 or 1.
    !       This means that the routing is limited to D8 (all the flux goes to one downstream cell),
    !        and partial refreezing is not allowed (i.e., btemp_weight_freeze = 1 everywhere).
    !       Thus, btemp_weight_freeze is not needed, but I kept it to keep the code similar to the subroutine above.
    !       We could make refreezing all-or-nothing (i.e., weights of either 0 or 1), but this leads to
    !        oscillations in bed temperature.
    !       I thought of rescaling flux_fraction and btemp_weight_freeze to largish i8 integers (e.g., 1000)
    !        to keep everything BFB, and then scaling back at the end. The problem is that this subroutine
    !        may need to be called repeatedly, and each scaling would lead to larger and larger integers
    !        that eventually exceed the i8 limit on integer size, ~10^(19).

    integer(i8), dimension(-1:1,-1:1,nx,ny), intent(in) ::  &
         flux_fraction           ! fraction of flux from a cell that flows downhill to each of 8 neighbors
                                 ! last two indices identify the source cell;
                                 ! 1st two indices give relative location of receiving cell

    integer(i8), dimension(nx,ny), intent(in) :: &
         btemp_weight_freeze     ! temperature-dependent weighting factor, favoring refreezing at frozen beds

    integer(i8), dimension(nx,ny), intent(inout) :: &
         bwatflx,              & ! on input: water flux (m^3/s * factor_bwatflx) to be routed to the margin or halo
                                 ! on output: flux routed to halo, to be routed further next time
         bwatflx_accum,        & ! cumulative bwatflx (m^3/s * factor_bwatflx) over multiple iterations
         bwatflx_refreeze_accum  ! cumulative bwatflx_refreeze (m^3/s * factor_bwatflx) over multiple iterations

    logical, intent(inout) :: &
         finished                ! initially F; set to T when all water has been routed as far as it can go

    ! Local variables

    integer :: i, j, k
    integer :: ii, jj, ip, jp

    ! Note: Some of the local variables are scaled by products of all three scale factors above.
    integer(i8), dimension(-1:1,-1:1,nx,ny)::  &
         bwatflx_halo            ! flux routed to halo cells
                                 ! last two indices identify the source cell;
                                 ! 1st two indices give relative location of receiving cell

    integer(i8), dimension(nx,ny) :: &
         bwatflx_refreeze,     & ! flux  saved for later refreezing; not routed further downstream
         sum_bwatflx_halo        ! bwatflx_halo summed over the first 2 indices

    integer(i8) :: &
         flx_thru,             & ! flux (m^3/s) that continues downstream
         global_halo_sum         ! global sum of water flux in halo cells

    real(dp), dimension(nx,ny)::  &
         bwatflx_dp, bwatflx_halo_dp, bwatflx_refreeze_dp  ! temporary dp versions of i8 arrays

    ! Initialize fluxes
    bwatflx_halo = 0
    bwatflx_refreeze = 0

    ! loop from high to low values on the local processor
    do k = nlocal, 1, -1

       ! Get i and j indices of current cell
       i = sorted_ij(k,1)
       j = sorted_ij(k,2)

       if (bwat_mask(i,j) == 1 .and. bwatflx(i,j) > 0.0d0) then

          ! Distribute the flux to downstream neighbors.
          ! Note: If btemp_weight_freeze = 1 everwhere, there is no refreezing.
          flx_thru = bwatflx(i,j) * btemp_weight_freeze(i,j)
          bwatflx_refreeze(i,j) = bwatflx(i,j) * (1 - btemp_weight_freeze(i,j))
          do jj = -1,1
             do ii = -1,1
                ip = i + ii
                jp = j + jj
                if (flux_fraction(ii,jj,i,j) > 0) then
                   if (halo_mask(ip,jp) == 1) then
                      bwatflx_halo(ii,jj,i,j) = flx_thru*flux_fraction(ii,jj,i,j)
                      if (verbose_bwat .and. this_rank==rtest .and. i==itest .and. j==jtest .and. count <= 2) then
                         write(iulog,*) 'Flux to halo, i, j, ii, jj, flux:', &
                              i, j, ii, jj, flx_thru*flux_fraction(ii,jj,i,j)
                      endif
                   elseif (local_mask(ip,jp) == 1) then
                      bwatflx(ip,jp) = bwatflx(ip,jp) + flx_thru*flux_fraction(ii,jj,i,j)
                      if (verbose_bwat .and. this_rank==rtest .and. i==itest .and. j==jtest .and. count <= 2) then
                         write(iulog,*) 'Flux to neighbor, i, j, ii, jj, flux:', &
                              i, j, ii, jj, flx_thru*flux_fraction(ii,jj,i,j)
                      endif
                   endif
                endif   ! flux_fraction > 0
             enddo  ! ii
          enddo  ! jj
       endif  ! bwat_mask = 1, bwatflx > 0
    enddo  ! loop from high to low

    ! Accumulate the fluxes in the output arrays
    bwatflx_accum = bwatflx_accum + bwatflx
    bwatflx_refreeze_accum = bwatflx_refreeze_accum + bwatflx_refreeze

    ! Compute the total bwatflx in halo cells
    do j = 1, ny
       do i = 1, nx
          sum_bwatflx_halo(i,j) = sum(bwatflx_halo(:,:,i,j))
       enddo
    enddo
    global_halo_sum = parallel_global_sum(sum_bwatflx_halo, parallel)

    ! If bwatflx_halo = 0 everywhere, then we are done.
    ! Where bwatflx_halo is nonzero, communicate it to the neighboring task.
    ! It will be routed further downstream the next time this subroutine is called.

    if (global_halo_sum > 0) then

       if (verbose_bwat .and. count <= 2) then
          if (this_rank == rtest) write(iulog,*) 'Before halo update, global_halo_sum (m^3/s):', global_halo_sum
       endif

       ! Reset bwatflx to zero for the halo transfer
       bwatflx = 0

       ! Communicate bmltflx_halo to the halo cells of neighboring processors
       call parallel_halo(bwatflx_halo(:,:,:,:), parallel)

       ! bmltflx_halo is now available in the halo cells of the local processor.
       ! Route downslope to the adjacent locally owned cells.
       ! These fluxes will be routed further downstream during the next iteration.
       do j = 2, ny-1
          do i = 2, nx-1
             if (halo_mask(i,j) == 1 .and. sum(bwatflx_halo(:,:,i,j)) > 0) then
                do jj = -1,1
                   do ii = -1,1
                      if (bwatflx_halo(ii,jj,i,j) > 0) then
                         ip = i + ii
                         jp = j + jj
                         if (local_mask(ip,jp) == 1) then
                            bwatflx(ip,jp) = bwatflx(ip,jp) + bwatflx_halo(ii,jj,i,j)
                            if (verbose_bwat .and. ip==itest .and. jp==jtest .and. this_rank==rtest .and. count <= 2) then
                               write(iulog,*) 'Nonzero bwatflx from halo, rank, i, j:', &
                                    this_rank, ip, jp, bwatflx_halo(ii,jj,i,j)
                            endif
                         endif
                      endif   !  bwatflx_halo > 0 to a local cell
                   enddo   ! ii
                enddo   ! jj
             endif   ! bwatflx_halo > 0 from this halo cell
          enddo   ! i
       enddo   ! j

    else

       finished = .true.   ! no water in halo cells to route further

    endif  ! global_halo_sum > 0

  end subroutine route_flux_to_margin_or_halo_integer8

!==============================================================

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !
  ! The following two subroutines perform an index-sort of an array.
  ! They are a GPL-licenced replacement for the Numerical Recipes routine indexx.
  ! They are not derived from any NR code, but are based on a quicksort routine by
  ! Michael Lamont (http://linux.wku.edu/~lamonml/kb.html), originally written
  ! in C, and issued under the GNU General Public License. Ian Rutt did the conversion
  ! to Fortran 90 and modified the algorithm to do an index sort.
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
