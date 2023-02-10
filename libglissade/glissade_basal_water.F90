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
   use glimmer_paramets, only: eps11, eps10, eps08
   use glimmer_physcon, only: rhoi, rhow, grav, scyr
   use glimmer_log
   use glide_types
   use cism_parallel, only: main_task, this_rank, nhalo, parallel_type, parallel_halo

   implicit none

   private
   public :: glissade_basal_water_init, glissade_calcbwat, glissade_bwat_flux_routing

!   logical, parameter :: verbose_bwat = .false.
   logical, parameter :: verbose_bwat = .true.

   integer, parameter :: pdiag = 3  ! range for diagnostic prints

contains

!==============================================================

  subroutine glissade_basal_water_init(model)

    use glimmer_paramets, only: thk0

    type(glide_global_type) :: model

    select case (model%options%which_ho_bwat)

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

        nx = size(bwat,1)
        ny = size(bwat,2)

        do j = 1, ny
           do i = 1, nx

              ! compute new bwat, given source (bmlt) and sink (drainage)
              ! Note: bmlt > 0 for ice melting. Freeze-on will reduce bwat.
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
       bmlt,                         &
       bwatflx,       head)

    ! This subroutine is a recoding of Jesse Johnson's steady-state water routing scheme in Glide.
    ! It has been parallelized for Glissade.
    !
    ! The subroutine returns the steady-state basal water flux, bwatflx,
    !  which reduces effective pressure N when which_ho_effecpress = HO_EFFECPRESS_BWATFLX.
    ! It should not be used with which_ho_effecpress = HO_EFFECPRESS_BWAT, since bwat is not returned.
    ! (See the comments below on bwat.)

    !TODO - Pass in a potential for basal freezing (bfrz_pot).
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
         bmlt                       ! basal melt rate (m/s)

    real(dp), intent(in) ::  &
         thklim                     ! minimum ice thickness for basal melt and hydropotential calculations (m)
                                    ! Note: This is typically model%geometry%thklim_temp

    integer, dimension(nx,ny), intent(in) ::  &
         bwat_mask,               & ! mask to identify cells through which basal water is routed;
                                    ! = 0 for floating and ocean cells; cells at global domain edge;
                                    !  and cells with thck = 0 and forced negative SMB
         floating_mask              ! = 1 if ice is present (thck > thklim) and floating, else = 0


    real(dp), dimension(nx,ny), intent(out) ::  &
         bwatflx,                 & ! basal water flux (m/yr)
         head                       ! hydraulic head (m)

    ! Local variables

    integer :: i, j, p

    real(dp), dimension(nx, ny) ::  &
         bwat,                    & ! diagnosed basal water depth (m), not used outside this module
         effecpress,              & ! effective pressure
         lakes                      ! difference between filled head and original head (m)

    ! parameters related to effective pressure
    real(dp), parameter :: &
         c_effective_pressure = 0.0d0    ! parameter in N = c/bwat; not currently used

    ! parameters related to subglacial fluxes
    ! The water flux q is given by Sommers et al. (2018), Eq. 5:
    !
    !           q = (b^3*g)/[(12*nu)*(1 + omega*Re)] * (-grad(h))
    !
    ! where q = basal water flux per unit width (m^2/s) = bwatflx/dx
    !       b = water depth (m) = bwat
    !       g = gravitational constant (m/s^2) = grad
    !      nu = kinematic viscosity of water (m^2/s)= visc_water
    !   omega = parameter controlling transition between laminar and turbulent flow
    !      Re = Reynolds number (large for turbulent flow)
    !       h = hydraulic head (m)
    !
    ! Note: In the equation above and the calculation below, bwatflx has units of m^3/s,
    !        i.e., volume per second entering and exiting a grid cell.
    !       For output, bwatflx has units of m/yr, i.e. volume per unit area per year entering and exiting a grid cell.
    !       With the latter convention, bwatflx is independent of grid resolution..
    !
    ! By default, we set Re = 0, which means the flow is purely laminar, as in Sommers et al. (2018), Eq. 6.

    ! Optionally, one or more of these parameters could be made a config parameter in the basal_hydro type
    real(dp), parameter ::  &
         visc_water = 1.787e-6,                       & ! kinematic viscosity of water (m^2/s); Sommers et al. (2018), Table 2
         omega_hydro = 1.0d-3,                        & ! omega (unitless) in Sommers et al (2018), Eq. 6
         Reynolds = 0.0d0                               ! Reynolds number (unitless), = 0 for pure laminar flow

    real(dp), parameter ::  &
         c_flux_to_depth = 1.0d0/((12.0d0*visc_water)*(1.0d0 + omega_hydro*Reynolds)), & ! proportionality coefficient in Eq. 6
         p_flux_to_depth = 2.0d0,                     & ! exponent for water depth; = 2 if q is proportional to b^3
         q_flux_to_depth = 1.0d0                        ! exponent for potential gradient; = 1 if q is linearly proportional to grad(h)

    ! WHL - debug fix_flats subroutine
    logical :: test_fix_flats = .false.
!!    logical :: test_fix_flats = .true.
    integer :: nx_test, ny_test
    real(dp), dimension(:,:), allocatable :: phi_test
    integer,  dimension(:,:), allocatable :: mask_test

    !WHL - debug
    !Note: This test works in serial, but does not work with parallel updates.
    !      To use it again, would need to comment out parallel calls in fix_flats.
    if (test_fix_flats) then

       ! Solve the example problem of Garbrecht & Martz (1997)
       ! Their problem is 7x7, but easier to solve if padded with a row of low numbers.

       nx_test = 9
       ny_test = 9
       allocate (phi_test(nx_test,ny_test))
       allocate (mask_test(nx_test,ny_test))

       mask_test = 1
       do j = 1, ny_test
          do i = 1, nx_test
             if (i == 1 .or. i == nx_test .or. j == 1 .or. j == ny_test) then
                mask_test(i,j) = 0
             endif
          enddo
       enddo

       phi_test(:,9) = (/ 1.d0, 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0, 1.d0 /)
       phi_test(:,8) = (/ 1.d0, 9.d0,9.d0,9.d0,9.d0,9.d0,9.d0,9.d0, 1.d0 /)
       phi_test(:,7) = (/ 1.d0, 9.d0,6.d0,6.d0,6.d0,6.d0,6.d0,9.d0, 1.d0 /)
       phi_test(:,6) = (/ 1.d0, 8.d0,6.d0,6.d0,6.d0,6.d0,6.d0,9.d0, 1.d0 /)
       phi_test(:,5) = (/ 1.d0, 8.d0,6.d0,6.d0,6.d0,6.d0,6.d0,9.d0, 1.d0 /)
       phi_test(:,4) = (/ 1.d0, 7.d0,6.d0,6.d0,6.d0,6.d0,6.d0,8.d0, 1.d0 /)
       phi_test(:,3) = (/ 1.d0, 7.d0,6.d0,6.d0,6.d0,6.d0,6.d0,8.d0, 1.d0 /)
       phi_test(:,2) = (/ 1.d0, 7.d0,7.d0,5.d0,7.d0,7.d0,8.d0,8.d0, 1.d0 /)
       phi_test(:,1) = (/ 1.d0, 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0, 1.d0 /)

       call fix_flats(&
            nx_test, ny_test,  &
            parallel,          &
            5,  5,   rtest,    &
            phi_test,          &
            mask_test)

       deallocate(phi_test, mask_test)

    endif

    if (verbose_bwat .and. this_rank == rtest) then
       print*, 'In glissade_bwat_flux_routing: rtest, itest, jtest =', rtest, itest, jtest
    endif

    ! Uncomment if the following fields are not already up to date in halo cells
!    call parallel_halo(thk,  parallel)
!    call parallel_halo(topg, parallel)
    call parallel_halo(bwat, parallel)
    call parallel_halo(bmlt, parallel)
    !TODO - Add bfrz?

    ! Compute effective pressure N.
    ! In the old Glimmer code, N was computed as a function of water depth by subroutine effective_pressure.
    ! Here, simply set N = 0 for the purpose of computing the hydraulic head.
    ! This approximation implies head = z_b + (rhoi/rhow) * H

!    call effective_pressure(&
!         bwat,                 &
!         c_effective_pressure, &
!         effecpress)

    effecpress(:,:) = 0.0d0

    ! Compute the hydraulic head

    call compute_head(&
         nx,     ny,    &
         thck,          &
         topg,          &
         effecpress,    &
         thklim,        &
         floating_mask, &
         head)

    p = pdiag

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'thck (m):'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i12)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f12.5)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'topg (m):'
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f12.5)',advance='no') topg(i,j)
          enddo
          write(6,*) ' '
       enddo
!       print*, ' '
!       print*, 'effecpress (Pa):'
!       write(6,*) ' '
!       do j = jtest+p, jtest-p, -1
!          write(6,'(i6)',advance='no') j
!          do i = itest-p, itest+p
!             write(6,'(f12.5)',advance='no') effecpress(i,j)
!          enddo
!          write(6,*) ' '
!       enddo
       print*, ' '
       print*, 'bmlt (m/yr):'
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f12.5)',advance='no') bmlt(i,j) * scyr
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'bwat_mask:'
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(i12)',advance='no') bwat_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Before fill: head (m):'
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f12.5)',advance='no') head(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Route basal water down the gradient of hydraulic head, giving a water flux
    ! TODO - Pass in bfrz_pot, return bfrz.

    call route_basal_water(&
         nx,      ny,            &
         dx,      dy,            &
         parallel,               &
         itest, jtest, rtest,    &
         flux_routing_scheme,    &
         head,                   &
         bmlt,                   &
         bwat_mask,              &
         bwatflx,                &
         lakes)

    call parallel_halo(bwatflx, parallel)

    ! Convert the water flux to a basal water depth
    ! Note: bwat is not passed out of this subroutine, for the following reason:
    ! In the thermal solve, the basal temperature is held at Tpmp wherever bwat > 0.
    ! This is appropriate when bwat is prognosed from local basal melting.
    ! For the flux-routing scheme, however, we can diagnose nonzero bwat beneath ice
    !  that is frozen to the bed (due to basal melting upstream).
    ! If passed to the thermal solver, this bwat can drive a sudden large increase in basal temperature.
    ! For this reason, the effective pressure is reduced based on bwatflx instead of bwat.
    ! If desired, we could pass out a diagnostic bwat field that would not affect basal temperature.

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
         bwat)

    ! Convert bwatflx units to m/yr for output
    bwatflx(:,:) = bwatflx(:,:) * scyr/(dx*dy)

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       write(6,*) 'Final bwatflx (m/yr):'
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f12.5)',advance='no') bwatflx(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Diagnosed bwat (mm):'
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f12.5)',advance='no') bwat(i,j) * 1000.d0
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine glissade_bwat_flux_routing

!==============================================================

  subroutine effective_pressure(&
       bwat,                    &
       c_effective_pressure,    &
       effecpress)

    ! Compute the effective pressure: the part of ice overburden not balanced by water pressure
    ! For now, this subroutine is not called; we just set effecpress = 0 for purposes of computing 'head'.
    ! TODO: Call calc_effecpress instead?

    real(dp),dimension(:,:),intent(in)  ::  bwat                  ! water depth
    real(dp)               ,intent(in)  ::  c_effective_pressure  ! constant of proportionality
    real(dp),dimension(:,:),intent(out) ::  effecpress            ! effective pressure

    ! Note: By default, c_effective_pressure = 0
    !       This implies N = 0; full support of the ice by water at the bed

    where (bwat > 0.d0)
        effecpress = c_effective_pressure / bwat
    elsewhere
        effecpress = 0.d0
    endwhere

  end subroutine effective_pressure

!==============================================================

  subroutine compute_head(&
       nx,      ny,   &
       thck,          &
       topg,          &
       effecpress,    &
       thklim,        &
       floating_mask, &
       head)

    !  Compute the hydraulic head as the bed elevation plus the scaled water pressure:
    !
    !     head = z_b + p_w / (rhow*g)
    !
    !  where z_b = bed elevation (m) = topg
    !        p_w = water pressure (Pa) = p_i - N
    !        p_i = ice overburden pressure = rhoi*g*H
    !          N = effective pressure (Pa) = part of overburden not supported by water
    !          H = ice thickness (m)
    !
    !  If we make the approximation p_w = p_i, then
    !
    !     head = z_b + (rhoi/rhow) * H

    implicit none

    ! Input/output variables

    integer, intent(in) ::  &
         nx, ny                   ! number of grid cells in each direction

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                  & ! ice thickness (m)
         topg,                  & ! bed elevation (m)
         effecpress               ! effective pressure (Pa)

    real(dp), intent(in) ::  &
         thklim                   ! minimum ice thickness for bmlt and head calculations

    integer, dimension(nx,ny), intent(in) ::  &
         floating_mask            ! = 1 if ice is present (thck > thklim) and floating, else = 0

    real(dp), dimension(nx,ny), intent(out) ::  &
         head                     ! hydraulic head (m)

    where (thck > thklim .and. floating_mask /= 1)
       head = topg + (rhoi/rhow)*thck - effecpress/(rhow*grav)
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
         bmlt,                   &
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
         bmlt                    ! basal melt beneath grounded ice (m/s)

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
         total_flux_in,  & ! total input flux (m^3/s), computed as sum of bmlt*dx*dy
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
          elseif (j > nhalo .and. j <= ny-nhalo .and. i > nhalo .and. i <= nx-nhalo) then
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

 ! Note: In an earlier code version, fix_flats was called here.
 !       It is not needed, however, if the fill_depressions scheme is run with epsilon > 0.

    ! Compute the lake depth
    lakes = head_filled - head

    ! Update head with the filled values
    head = head_filled

    p = pdiag
    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'After fill: head (m):'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i12)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f12.5)',advance='no') head(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'lakes (m):'
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f12.5)',advance='no') lakes(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Sort heights.
    ! The sorted_ij array stores the i and j index for each locally owned cell, from lowest to highest value.

    call sort_heights(&
         nx,    ny,    nlocal,  &
         itest, jtest, rtest,   &
         head,  sorted_ij)

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
    ! Note: This subroutine conserves water only if bmlt >= 0 everywhere.
    !       One way to account for refreezing would be to do the thermal calculation after
    !        computing bwat in this subroutine.  At that point, refreezing would take away
    !        from the bwat computed here.  In the next time step, positive values of bmlt
    !        would provide a new source for bwat.
    ! In other words, the sequence would be:
    ! (1) Ice transport and calving
    ! (2) Basal water routing: apply bmlt and diagnose bwat
    ! (3) Vertical heat flow:
    !     (a) compute bmlt
    !     (b) use bmlt < 0 to reduce bwat
    !     (c) save bmlt > 0 for the next time step (and write to restart)
    ! (4) Diagnose velocity

    bwatflx = 0.0d0
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          bwatflx(i,j) = bmlt(i,j) * dx * dy
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
          print*, ' '
          print*, 'Current bwatflx_accum, rank, i, j =', this_rank, itest, jtest
          write(6,*) ' '
          do j = jtest+p, jtest-p, -1
             write(6,'(i6)',advance='no') j
             do i = itest-p, itest+p
                write(6,'(e12.3)',advance='no') bwatflx_accum(i,j)
             enddo
             write(6,*) ' '
          enddo
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
!                print*, 'Nonzero bwatflx_halo, count, rank, i, j, sum_bwatflx_halo:', &
!                     count, this_rank, i, j, sum_bwatflx_halo(i,j)
!                call parallel_globalindex(i, j, iglobal, jglobal, parallel)
!                print*, '     iglobal, jglobal:', iglobal, jglobal
!             endif
          enddo
       enddo
       global_flux_sum = parallel_global_sum(sum_bwatflx_halo, parallel)

       if (verbose_bwat .and. this_rank == rtest .and. count <= 2) then
          print*, ' '
          print*, 'Before halo update, sum of bwatflx_halo:', global_flux_sum
          print*, ' '
          print*, 'sum_bwatflx_halo:'
          write(6,*) ' '
          do j = jtest+p, jtest-p, -1
             write(6,'(i6)',advance='no') j
             do i = itest-p, itest+p
                write(6,'(e12.3)',advance='no') sum_bwatflx_halo(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'rank, i, j, bwatflx_halo:'
          i = itest
          j = jtest
          write(6, '(3i5,9e10.3)') this_rank, i, j, bwatflx_halo(:,:,i,j)
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
       print*, ' '
       print*, 'Total output basal melt flux (m^3/s):', total_flux_out
       print*, 'Difference between input and output =', total_flux_in - total_flux_out
    endif

    ! Found that an error threshold of eps11 is not quite large enough for Antarctic 8km runs.
    if (total_flux_in > 0.0d0) then
       err = abs(total_flux_in - total_flux_out)
       if (err > eps10) then
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
         bwat)

    !  Assuming that the flow is steady state, this function simply solves
    !               flux = depth * velocity
    !  for the depth, assuming that the velocity is a function of depth,
    !  and pressure potential. This amounts to assuming a Weertman film,
    !  or Manning flow, both of which take the form of a constant times water
    !  depth to a power, times grad(head) to a power.

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
         bwat_mask               ! mask to identify cells through which basal water is routed;
                                 ! = 1 where ice is present and not floating

    real(dp), dimension(nx,ny), intent(out)::  &
         bwat                     ! water depth

    ! Local variables

    integer :: i, j, p

    real(dp), dimension(nx,ny) ::  &
         grad_head                ! gradient of hydraulic head (m/m), averaged to cell centers

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
    p = pdiag
    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'grad_head:'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i12)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f12.5)',advance='no') grad_head(i,j)
          enddo
          write(6,*) ' '
       enddo
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
    ! Jesse's Glimmer code has this:
    !       bwat = ( bwatflx / (c_flux_to_depth * scyr * dy * grad_wphi**q_flux_to_depth) ) ** exponent
    ! which is missing the grav term and seems to have an extra scyr term.
    ! Also, c_flux_to_depth = 1 / (12 * 1.6d-3) in Jesse's code.  Note exponent of d-3 instead of d-6 for nu.
    !
    ! Note: Assuming dx = dy
    ! TODO: Modify for the case dx /= dy?

    where (grad_head /= 0.d0 .and. bwat_mask == 1)
       bwat = ( bwatflx / (c_flux_to_depth * grav * dy * grad_head**q_flux_to_depth) ) ** p_exponent
    elsewhere
       bwat = 0.d0
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
    ! (The boundary consists of cells with phi_mask = 0, adjacent to cells with phi_mask = 1.)
    ! Loop through the domain.  For each cell c, with value phi(c) not yet fixed as a known value,
    !  compute phi_min8(n), the current minimum of phi in the 8 neighbor cells.
    ! If phi_in(c) > phi_min8(n) + eps, then set phi(c) = phi_in(c) and mark that cell as having a known value,
    !  since phi(c) cannot go any lower.  Here, eps is a small number greater than zero.
    ! If phi_in(c) < phi_min8(n) + eps, but phi(c) > phi_min8(c) + eps, set phi(c) = phi_min8(n) + eps.
    !  Do not mark the cell as having a known value, because it might be lowered further.
    ! Continue until no further lowering of phi is possible.  At that point, phi = phi_out.

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
         epsilon_edge = 1.d-4,          & ! small increment in phi to avoid flat regions, applied to edge neighbors
         epsilon_diag = 1.d-4*sqrt(2.d0)  ! small increment in phi to avoid flat regions, applied to diagonal neighbors

    !WHL - Typically, it takes ~10 iterations to fill all depressions on a large domain.
    integer, parameter :: count_max = 100

    logical, parameter :: verbose_depression = .false.
!!    logical, parameter :: verbose_depression = .true.

    ! Initial halo update, in case phi_in is not up to date in halo cells
    call parallel_halo(phi_in, parallel)

    ! Initialize phi to a large value
    where (phi_mask == 1)
       phi = big_number
       known = .false.
    elsewhere
       phi = 0.0d0
       known = .true.
    endwhere

    ! Set phi = phi_in for boundary cells.
    ! A boundary cell is a cell with phi_mask = 0, adjacent to one or more cells with phi_mask = 1.

    do j = 2, ny-1
       do i = 2, nx-1
          if (phi_mask(i,j) == 0) then
             if (phi_mask(i-1,j+1)==1 .or. phi_mask(i,j+1)==1 .or. phi_mask(i+1,j+1)==1 .or. &
                 phi_mask(i-1,j)  ==1               .or.           phi_mask(i+1,j)  ==1 .or. &
                 phi_mask(i-1,j-1)==1 .or. phi_mask(i,j-1)==1 .or. phi_mask(i+1,j-1)==1) then
                phi(i,j) = phi_in(i,j)
                known(i,j) = .true.
             endif
          endif
       enddo
    enddo

    ! The resulting field phi applies to locally owned cells and one layer of halo cells.
    ! A halo update brings it up to date in all halo cells.

    call parallel_halo(phi, parallel)

    p = pdiag

    if (verbose_depression .and. this_rank == rtest) then
       print*, ' '
       print*, 'Initial phi:'
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(e11.4)',advance='no') phi(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    count = 0
    finished = .false.

    do while (.not.finished)

       count = count + 1
       local_lowered = 0

       if (verbose_depression .and. this_rank == rtest) then
          write(6,*) ' '
          print*, 'fill_depressions, count =', count
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

       if (verbose_depression .and. this_rank == rtest) then
          print*, ' '
          print*, 'New phi:'
          do j = jtest+p, jtest-p, -1
             write(6,'(i6)',advance='no') j
             do i = itest-p, itest+p
                write(6,'(f11.4)',advance='no') phi(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       ! If one or more cells was lowered, then repeat; else exit the local loop.

       global_lowered = parallel_reduce_sum(local_lowered)

       if (global_lowered == 0) then
          finished = .true.
          if (verbose_depression .and. this_rank == rtest) then
             print*, 'finished lowering'
          endif
       else
          finished = .false.
          if (verbose_depression .and. this_rank == rtest) then
             print*, 'cells lowered on this iteration:', global_lowered
          endif
          call parallel_halo(phi, parallel)
       endif

       if (count > count_max) then
          call write_log('Hydrology error, exceeded max number of global iterations', GM_FATAL)
       endif

    enddo  ! finished

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'Filled depressions, count =', count
    endif

  end subroutine fill_depressions

!==============================================================

  !TODO - Delete the fix_flats subroutine; replaced by fill_depressions

  subroutine fix_flats(&
       nx,    ny,            &
       parallel,             &
       itest, jtest, rtest,  &
       phi,                  &
       phi_mask)

    ! Add a small increment to flat regions in the input field phi,
    !  so that all cells have a downhill outlet.
    !
    ! Use the algorithm of Garbrecht & Mertz:
    ! Garbrecht, J., and L. W. Mertz (1997), The assignment of drainage direction
    !    over flat surfaces in raster digital elevation models, J. Hydrol., 193,
    !    204-213.
    !
    ! Note: This subroutine is not currently called, because the depression-filling algorithm
    !        above is configured not to leave any flats.
    !       I am leaving it here in case it is useful for debugging.

    use cism_parallel, only: parallel_global_sum

    implicit none

    ! Input/output variables

    integer, intent(in) ::  &
         nx, ny,                & ! number of grid cells in each direction
         itest, jtest, rtest      ! coordinates of diagnostic point

    type(parallel_type), intent(in) ::  &
         parallel                 ! info for parallel communication

    real(dp), dimension(nx,ny), intent(inout) :: &
         phi                      ! input field with flat regions to be fixed

    integer, dimension(nx,ny), intent(in) ::  &
         phi_mask                 ! = 1 where any flat regions of phi will need to be fixed, else = 0
                                  ! corresponds to the grounded ice sheet (bmlt_mask) for the flux-routing problem

    ! Local variables --------------------------------------

    real(dp), dimension(nx,ny) ::  &
         phi_input,             & ! input value of phi, before any corrections
         phi_new,               & ! new value of phi, after incremental corrections
         dphi1,                 & ! sum of increments applied in step 1
         dphi2                    ! sum of increments applied in step 2

    integer, dimension(nx,ny) :: &
         flat_mask,             & ! = 1 for cells with phi_mask = 1 and without a downslope gradient, else = 0
         flat_mask_input,       & ! flat_mask as computed from phi_input
         n_uphill,              & ! number of uphill neighbors for each cell, as computed from input phi
         n_downhill,            & ! number of downhill neighbors for each cell, as computed from input phi
         incremented_mask,      & ! = 1 for cells that have already been incremented (in step 2)
         unincremented_mask,    & ! = 1 for cells in input flat regions, not yet incremented
         incremented_neighbor_mask  ! = 1 for cells that have not been incremented, but have an incremented neighbor

    integer :: &
         global_sum               ! global sum of cells meeting a mask criterion

    logical :: finished           ! true when an iterative loop has finished

    real(dp), parameter :: &
         phi_increment = 2.0d-5   ! fractional increment in phi (Garbrecht & Martz use 2.0e-5)

    integer :: i, j, ii, jj, ip, jp, p
    integer :: count
    integer, parameter :: count_max = 100

    ! Uncomment if the input fields are not up to date in halos
!    call parallel_halo(phi, parallel)
!    call parallel_halo(phi_mask, parallel)

    p = pdiag

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'In fix_flats, rtest, itest, jtest =', rtest, itest, jtest
       print*, ' '
       print*, 'input phi:'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i12)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f12.5)',advance='no') phi(i,j)
          enddo
          write(6,*) ' '
       enddo
       write(6,*) ' '
       print*, 'phi_mask:'
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(i12)',advance='no') phi_mask(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! initialize

    phi_input = phi

    ! For use in Step 2, count the uphill and downhill neighbors of each cell.

    n_uphill = 0
    n_downhill = 0

    do j = 2, ny-1
       do i = 2, nx-1
          if (phi_mask(i,j) == 1) then
             do jj = -1,1
                do ii = -1,1
                   ! If this is the centre point, ignore
                   if (ii == 0 .and. jj == 0) then
                      continue
                   else  ! check for nonzero elevation gradients
                      ip = i + ii
                      jp = j + jj
                      if (phi(ip,jp) - phi(i,j) > eps11) then  ! uphill neighbor
                         n_uphill(i,j) = n_uphill(i,j) + 1
                      elseif (phi(i,j) - phi(ip,jp) > eps11) then  ! downhill neighbor
                         n_downhill(i,j) = n_downhill(i,j) + 1
                      endif
                   endif
                enddo   ! ii
             enddo   ! jj
          endif   ! phi_mask = 1
       enddo   ! i
    enddo   ! j

    ! Identify the flat regions in the input field.
    ! This includes all cells with phi_mask = 1 and without downslope neighbors.
    ! The resulting flat_mask is valid in all cells except the outer halo.

    call find_flats(&
         nx,    ny,           &
         itest, jtest, rtest, &
         phi_input,           &
         phi_mask,            &
         flat_mask_input)

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'n_uphill:'
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(i10)',advance='no') n_uphill(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'n_downhill:'
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(i10)',advance='no') n_downhill(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'input flat_mask:'
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(i10)',advance='no') flat_mask_input(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Step 1: Gradient toward lower terrain

    dphi1 = 0.0d0
    flat_mask = flat_mask_input
    finished = .false.
    count = 0

    ! Increment phi in all cells with flat_mask = 1 (no downslope gradient).
    ! Repeat until all cells on the global grid have a downslope gradient.

    do while(.not.finished)

       count = count + 1
       if (verbose_bwat .and. this_rank == rtest) then
          print*, ' '
          print*, 'step 1, count =', count
       endif

       where (flat_mask == 1)
          dphi1 = dphi1 + phi_increment
       endwhere

       call parallel_halo(dphi1, parallel)

       phi_new = phi_input + dphi1

       ! From the original flat region, identify cells that still have no downslope gradient.
       ! The resulting flat_mask is valid in all cells except the outer halo.

       call find_flats(&
            nx,    ny,             &
            itest, jtest, rtest,   &
            phi_new,               &
            flat_mask_input,       &
            flat_mask)

       if (verbose_bwat .and. this_rank == rtest) then
          print*, ' '
          print*, 'Updated dphi1/phi_increment:'
          do j = jtest+p, jtest-p, -1
             write(6,'(i6)',advance='no') j
             do i = itest-p, itest+p
                write(6,'(f10.1)',advance='no') dphi1(i,j)/ phi_increment
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'Updated flat_mask:'
          do j = jtest+p, jtest-p, -1
             write(6,'(i6)',advance='no') j
             do i = itest-p, itest+p
                write(6,'(i10)',advance='no') flat_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       ! Compute the number of cells in the remaining flat regions on the global grid.
       ! If there are no such cells, then exit the loop.

       global_sum = parallel_global_sum(flat_mask, parallel)

       if (verbose_bwat .and. this_rank == rtest) then
          print*, 'global sum of flat_mask =', global_sum
       endif

       if (global_sum > 0) then
          finished = .false.
       else
          finished = .true.
       endif

       if (count > count_max) then
          call write_log('Hydrology error: abort in step 1 of fix_flats', GM_FATAL)
       endif

    enddo   ! step 1 finished

    ! Step 2: Gradient away from higher terrain

    dphi2 = 0.0d0
    incremented_mask = 0
    finished = .false.
    count = 0

    ! In the first pass, increment the elevation in all cells of the input flat region that are
    !  adjacent to higher terrain and have no adjacent downhill cell.
    ! The resulting dphi2 and incremented_mask are valid in all cells except the outer halo.
    ! Need a halo update for incremented_mask to compute incremented_neighbor_mask below.

    do j = 2, ny-1
       do i = 2, nx-1
          if (flat_mask_input(i,j) == 1) then
             if (n_uphill(i,j) >= 1 .and. n_downhill(i,j) == 0) then
                dphi2(i,j) = dphi2(i,j) + phi_increment
                incremented_mask(i,j) = 1
             endif
          endif
       enddo
    enddo

    call parallel_halo(dphi2, parallel)
    call parallel_halo(incremented_mask, parallel)

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'step 2, input flat_mask:'
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(i10)',advance='no') flat_mask_input(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Updated dphi2/phi_increment'
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.1)',advance='no') dphi2(i,j)/phi_increment
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Compute the number of cells incremented in the first pass.
    ! If no cells are incremented, then skip step 2.
    ! This will be the case if the flat region lies at the highest elevation in the domain.

    global_sum = parallel_global_sum(incremented_mask, parallel)

    if (global_sum == 0) then
       if (verbose_bwat .and. this_rank == rtest) then
          print*, ' '
          print*, 'No cells to increment; skip step 2'
       endif
       finished = .true.
    endif

    ! In subsequent passes, increment the elevation in the following cells:
    !  (1) all cells that have been previously incremented; and
    !  (2) all cells in the input flat region that have not been previously incremented,
    !      are adjacent to an incremented cell, and are not adjacent to a cell downhill
    !      from the input flat region.
    ! Repeat until no unincremented cells remain on the input flat region.
    ! Note: This iterated loop uses flat_mask_input, which is not incremented.

    do while(.not.finished)

       count = count + 1
       if (verbose_bwat .and. this_rank == rtest) then
          print*, ' '
          print*, 'step 2, count =', count
       endif

       ! Identify cells that have not been incremented, but are adjacent to incremented cells
       ! The resulting incremented_neighbor mask is valid in all cells except the outer halo.

       incremented_neighbor_mask = 0
       do j = 2, ny-1
          do i = 2, nx-1
             if (flat_mask_input(i,j) == 1 .and. incremented_mask(i,j) == 0) then
                do jj = -1,1
                   do ii = -1,1
                      ! If this is the centre point, ignore
                      if (ii == 0 .and. jj == 0) then
                         continue
                      else  ! check for an incremented neighbor
                         ip = i + ii
                         jp = j + jj
                         if (incremented_mask(ip,jp) == 1) then
                            incremented_neighbor_mask(i,j) = 1
                         endif
                      endif
                   enddo   ! ii
                enddo   ! jj
             endif   ! flat_mask = 1 and incremented = F
          enddo   ! i
       enddo   ! j

       ! Increment cells of type (1) and (2)
       ! Note: n_downhill was computed before step 1.

       do j = 2, ny-1
          do i = 2, nx-1
             if (incremented_mask(i,j) == 1) then
                dphi2(i,j) = dphi2(i,j) + phi_increment
             elseif (n_downhill(i,j) == 0 .and. incremented_neighbor_mask(i,j) == 1) then
                dphi2(i,j) = dphi2(i,j) + phi_increment
                incremented_mask(i,j) = 1
             endif
          enddo
       enddo

       call parallel_halo(dphi2, parallel)
       call parallel_halo(incremented_mask, parallel)

       if (verbose_bwat .and. this_rank == rtest) then
          print*, ' '
          print*, 'incremented_neighbor_mask:'
          do j = jtest+p, jtest-p, -1
             write(6,'(i6)',advance='no') j
             do i = itest-p, itest+p
                write(6,'(i10)',advance='no') incremented_neighbor_mask(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, 'Updated dphi2/phi_increment'
          do j = jtest+p, jtest-p, -1
             write(6,'(i6)',advance='no') j
             do i = itest-p, itest+p
                write(6,'(f10.1)',advance='no') dphi2(i,j)/phi_increment
             enddo
             write(6,*) ' '
          enddo
       endif

       ! Compute the number of cells in the input flat region that have not been incremented.
       ! If all the flat cells have been incremented, then exit the loop.

       where (flat_mask_input == 1 .and. incremented_mask == 0)
          unincremented_mask = 1
       elsewhere
          unincremented_mask = 0
       endwhere
       global_sum = parallel_global_sum(unincremented_mask, parallel)


       if (global_sum > 0) then
          if (verbose_bwat .and. this_rank == rtest) then
             print*, 'number of flat cells not yet incremented =', global_sum
          endif
          finished = .false.
       else
          finished = .true.
       endif

       if (count > count_max) then
          call write_log('Hydrology error: abort in step 2 of fix_flats', GM_FATAL)
       endif

    enddo   ! step 2 finished


    ! Step 3:

    ! Add the increments from steps 1 and 2
    ! The result is a surface with gradients both toward lower terrain and away from higher terrain.
    ! No halo update is needed here, since dphi1 and dphi2 have been updated in halos.

    phi = phi_input + dphi1 + dphi2

    ! Check for cells with flat_mask = 1 (no downslope gradient).
    ! Such cells are possible because of cancelling dphi1 and dphi2.

    count = 0
    finished = .false.

    do while (.not.finished)

       count = count + 1
       if (verbose_bwat .and. this_rank == rtest) then
          print*, ' '
          print*, 'step 3, count =', count
       endif

       ! Identify cells without downslope neighbors.
       ! The resulting flat_mask is valid in all cells except the outer halo.

       call find_flats(&
            nx,    ny,           &
            itest, jtest, rtest, &
            phi,                 &
            phi_mask,            &
            flat_mask)

       ! Add a half increment to any cells without downslope neighbors
       where (flat_mask == 1)
          phi = phi + 0.5d0 * phi_increment
       endwhere

       call parallel_halo(phi, parallel)

       ! Compute the number of cells without downslope neighbors.
       ! If there are no such cells, then exit the loop.

       global_sum = parallel_global_sum(flat_mask, parallel)

       if (verbose_bwat .and. this_rank == rtest) then
          print*, 'global sum of flat_mask =', global_sum
       endif

       if (global_sum > 0) then
          finished = .false.
       else
          finished = .true.
       endif

       if (count > count_max) then
          call write_log('Hydrology error: abort in step 3 of fix_flats', GM_FATAL)
       endif

    enddo   ! step 3 finished

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'Final phi:'
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.3)',advance='no') phi(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine fix_flats

!==============================================================

  subroutine find_flats(&
       nx,     ny,             &
       itest,  jtest,  rtest,  &
       phi,    phi_mask,       &
       flat_mask)

    ! Compute a mask that = 1 for cells in flat regions.
    ! These are defined as cells with phi_mask = 1 and without a downslope gradient.
    ! Note: This subroutine is not currently called.  See the comment in fix_flats.

    ! Input/output arguments

    integer, intent(in) ::  &
         nx, ny,               & ! number of grid cells in each direction
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), dimension(nx,ny), intent(inout) :: &
         phi                     ! elevation field with potential flat regions

    integer, dimension(nx,ny), intent(in) :: &
         phi_mask               ! = 1 for cells in the region where flats need to be identified

    integer, dimension(nx,ny), intent(out) :: &
         flat_mask               ! = 1 for cells with phi_mask = 1 and without a downslope gradient

    ! Local variables

    integer :: i, j, ii, jj, ip, jp

    where (phi_mask == 1)
       flat_mask = 1   ! assume flat until shown otherwise
    elsewhere
       flat_mask = 0
    endwhere

    ! Identify cells that have no downslope neighbors, and mark them as flat.

    do j = 2, ny-1
       do i = 2, nx-1
          if (phi_mask(i,j) == 1) then
             !TODO - Add an exit statement?
             do jj = -1,1
                do ii = -1,1
                   ! If this is the centre point, ignore
                   if (ii == 0 .and. jj == 0) then
                      continue
                   else  ! check for a downslope gradient
                      ip = i + ii
                      jp = j + jj
                      if (phi(i,j) - phi(ip,jp) > eps11) then
                         flat_mask(i,j) = 0
                      endif
                   endif
                enddo   ! ii
             enddo   ! jj
          endif   ! phi_mask = 1
       enddo   ! i
    enddo   ! j

!    call parallel_halo(flat_mask, parallel)

  end subroutine find_flats

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

          endif   ! flux_routing_scheme: D8, Dinf, FD8

          if (verbose_bwat .and. this_rank == rtest .and. i == itest .and. j == jtest) then
             print*, 'sum_slope:', sum_slope
             print*, 'slope(:, 1):', slope(:, 1)
             print*, 'slope(:, 0):', slope(:, 0)
             print*, 'slope(:,-1):', slope(:,-1)
             print*, 'flux_fraction(:, 1,i,j):', flux_fraction(:, 1,i,j)
             print*, 'flux_fraction(:, 0,i,j):', flux_fraction(:, 0,i,j)
             print*, 'flux_fraction(:,-1,i,j):', flux_fraction(:,-1,i,j)
          endif

          ! bug check - make sure fractions add to 1
          sum_frac = 0.0d0
          do jj = -1,1
             do ii = -1,1
                sum_frac = sum_frac + flux_fraction(ii,jj,i,j)
             enddo
          enddo
          if (abs(sum_frac - 1.0d0) > eps11) then
             ! Do a fatal abort?
             print*, 'ERROR: sum_frac /= 1: r, i, j, sum:', this_rank, i, j, sum_frac
          endif

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
