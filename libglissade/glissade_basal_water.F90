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

!TODO - Test and parallelize Jesse's water-routing code.
!       Currently supported only for serial Glide runs, in module glide_bwater.F90

module glissade_basal_water

   use glimmer_global, only: dp
   use glimmer_paramets, only: eps11
   use glimmer_log
   use glide_types
   use parallel_mod, only: main_task, this_rank, parallel_type, parallel_halo

   implicit none

   private
   public :: glissade_basal_water_init, glissade_calcbwat, glissade_bwat_flux_routing

   logical, parameter :: verbose_bwat = .true.

   integer, parameter :: pdiag = 5  ! range for diagnostic prints

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

    use glimmer_physcon, only: rhow, scyr
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
       nx,            ny,      &
       dx,            dy,      &
       itest, jtest,  rtest,   &
       flux_routing_scheme,    &
       thklim,                 &
       thck,                   &
       topg,                   &
       bmlt,                   &
       floating_mask,          &
       bwat,                   &
       bwatflx,                &
       head)

    ! This subroutine is a recoding of Jesse Johnson's steady-state water routing scheme in Glide.
    ! Needs to be parallelized for Glissade.

    use glimmer_physcon, only: scyr
    use glimmer_log
    use parallel_mod, only: tasks   ! while code is serial only

    ! Input/output arguments

    integer, intent(in) ::  &
         nx, ny,             &      ! number of grid cells in each direction
         itest, jtest, rtest        ! coordinates of diagnostic point

    real(dp), intent(in) ::  &
         dx, dy                     ! grid cell size (m)

    integer, intent(in) ::  &
         flux_routing_scheme        ! flux routing scheme: D8, Dinf or FD8; see subroutine route_basal_water

    real(dp), intent(in) ::  &
         thklim                     ! minimum ice thickness for basal melt and hydropotential calculations (m)

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,               &      ! ice thickness (m)
         topg,               &      ! bed topography (m)
         bmlt                       ! basal melt rate (m/s)

    integer, dimension(nx,ny), intent(in) :: &
         floating_mask              ! = 1 if ice is present (thck > thklim) and floating, else = 0

    real(dp), dimension(nx,ny), intent(inout) ::  &
         bwat                       ! basal water depth (m)

    real(dp), dimension(nx,ny), intent(out) ::  &
         bwatflx,                 & ! basal water flux (m^3/s)
         head                       ! hydraulic head (m)

    ! Local variables

    integer :: i, j, p

    !TODO - Make effecpress in/out?
    real(dp), dimension(nx, ny) ::  &
         effecpress,    & ! effective pressure
         lakes            ! difference between filled head and original head (m)

    integer, dimension(nx,ny) :: &
         bwat_mask        ! mask to identify cells through which basal water is routed;
                          ! = 1 if ice is present (thck > thklim) and not floating, else = 0

    ! parameters related to effective pressure
    real(dp), parameter :: &
         c_effective_pressure = 0.0d0                   ! for now estimated as N = c/bwat

    ! parameters related to subglacial fluxes
    ! The basal water flux is given by Sommers et al. (2018), Eq. 5:
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
            5,  5,   rtest,    &
            phi_test,          &
            mask_test)

       deallocate(phi_test, mask_test)

    endif

    !WHL - debug
    if (main_task) print*, 'In glissade_bwat_flux_routing: rtest, itest, jtest =', rtest, itest, jtest

    if (tasks > 1) then
       call write_log('Flux routing not yet supported for tasks > 1', GM_FATAL)
    endif


    ! Compute effective pressure N as a function of water depth
    call effective_pressure(&
         bwat,                 &
         c_effective_pressure, &
         effecpress)

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
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.3)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'topg (m):'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.3)',advance='no') topg(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'effecpress (Pa):'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.3)',advance='no') effecpress(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'bmlt (m/yr):'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.3)',advance='no') bmlt(i,j) * scyr
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Before fill: head (m):'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.3)',advance='no') head(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Compute a mask: = 1 where ice is present and not floating
    where (thck > thklim .and. floating_mask == 0)
       bwat_mask = 1
    elsewhere
       bwat_mask = 0
    endwhere

    ! Route basal water down the gradient of hydraulic head, giving a water flux
    call route_basal_water(&
         nx,      ny,            &
         dx,      dy,            &
         itest, jtest, rtest,    &
         flux_routing_scheme,    &
         head,                   &
         bmlt,                   &
         bwat_mask,              &
         bwatflx,                &
         lakes)

    ! Convert the water flux to a basal water depth
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

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'bwatflx (m^3/s):'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.3)',advance='no') bwatflx(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'bwat (mm):'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.3)',advance='no') bwat(i,j) * 1000.d0
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
    ! TODO: Try c_effective_pressure > 0, or call calc_effecpress instead

    real(dp),dimension(:,:),intent(in)  ::  bwat                  ! water depth
    real(dp)               ,intent(in)  ::  c_effective_pressure  ! constant of proportionality
    real(dp),dimension(:,:),intent(out) ::  effecpress            ! effective pressure

    ! Note: By default, c_effective_pressure = 0
    !       This implies N = 0; full support of the ice by water at the bed
    !       Alternatively, could call the standard glissade subroutine, calc_effective_pressure

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
    !  If we make the approximation p_w =~ p_i, then
    !
    !     head =~ z_b + (rhoi/rhow) * H

    use glimmer_physcon, only : rhoi, rhow, grav
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
         itest, jtest, rtest,    &
         flux_routing_scheme,    &
         head,                   &
         bmlt,                   &
         bwat_mask,              &
         bwatflx,                &
         lakes)

    ! Route water from the basal melt field to its destination, recording the water flux along the route.
    ! Water flow direction is determined according to the gradient of the hydraulic head.
    ! For the algorithm to function properly, surface depressions must be filled,
    !  so that all cells have an outlet to the ice sheet margin.
    !> This results in the lakes field, which is the difference between the filled head and the original head.
    !> The method used is by Quinn et. al. (1991).
    !>
    !> Based on code by Jesse Johnson (2005), adapted from the glimmer_routing file by Ian Rutt.

    !TODO: This is a serial subroutine.
    !      To run in Glissade, we need to add a global gather/scatter.
    !      Ultimately, the goal is to make it fully parallel.

    implicit none

    integer, intent(in) ::  &
         nx, ny,               & ! number of grid cells in each direction
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), intent(in) ::  &
         dx, dy                  ! grid spacing in each direction (m)

    integer, intent(in) ::  &
         flux_routing_scheme     ! flux routing scheme: D8, Dinf or FD8
                                 ! D8: Flow is downhill toward the single cell with the lowest elevation.
                                 ! Dinf: Flow is downhill toward the two cells with the lowest elevations.
                                 ! FD8: Flow is downhill toward all cells with lower elevation.
                                 ! D8 scheme gives the narrowest flow, and FD8 gives the most diffuse flow.

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

    integer :: i, j, k, nn, ii, jj, ip, jp
    integer :: i1, j1, i2, j2, itmp, jtmp
    integer :: p

    integer,  dimension(:,:), allocatable ::  &
         sorted            ! i and j indices of all cells, sorted from low to high potential

    real(dp), dimension(nx,ny) ::  &
         head_filled       ! head after depressions are filled (m)

    integer, dimension(nx,ny) ::  &
         flats             !

    real(dp), dimension(-1:1,-1:1) ::  &
         dists,         &  ! distance (m) to adjacent grid cell
         slope             ! slope of head between adjacent grid cells

    real(dp) ::  &
         slope1,        &  ! largest value of slope array
         slope2,        &  ! second largest value of slope array
         sum_slope,     &  ! slope1 + slope2
         slope_tmp         ! temporary slope value

    logical :: flag

    real(dp) :: &
         total_flux_in, &  ! total input flux (m^3/s), computed as sum of bmlt*dx*dy
         total_flux_out, & ! total output flux (m^3/s), computed as sum of bwatflx at ice margin
         flux_unrouted     ! total flux (m^3/s) that is not routed downhill (should = 0)

    integer, dimension(nx,ny) ::  &
         margin_mask       ! = 1 for cells at the margin, as defined by bwat_mask


    ! Compute distances to adjacent grid cells for slope determination

    dists(-1,:) = (/ sqrt(dx**2 + dy**2), dy, sqrt(dx**2 + dy**2) /)
    dists(0,:) = (/ dx, 0.0d0, dx /)
    dists(1,:) = dists(-1,:)

    ! Allocate local arrays

    nn = nx*ny   ! For parallel code, change to locally owned cells only
    allocate(sorted(nn,2))

    ! Initialize the filled field
    head_filled = head

    ! Fill depressions in head, so that no interior cells are sinks
    call fill_depressions(&
         nx,          ny,    &
         head_filled,        &
         bwat_mask)


    ! Raise the head slightly in flat regions, so that all cells have downslope outlets

    call fix_flats(&
         nx,    ny,            &
         itest, jtest, rtest,  &
         head_filled,          &
         bwat_mask)

    lakes = head_filled - head

    p = pdiag
    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'After fill: head_filled (m):'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.3)',advance='no') head_filled(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'lakes (m):'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.3)',advance='no') lakes(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    ! Update head with the filled values
    head = head_filled

    ! Sort heights.
    ! The 'sorted' array contains the i and j index for each cell, from lowest to highest value of the filled potential.
    call heights_sort(&
         nx,           ny,     &
         itest, jtest, rtest,  &
         head,         sorted)

    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'sorted, from the top:'
       do k = nx*ny, nx*ny-10, -1
          i = sorted(k,1)
          j = sorted(k,2)
          print*, i, j, head(i,j)
       enddo
    endif

    ! Initialise the water flux with the local basal melt, which will then be redistributed.
    ! Multiply by area, so units are m^3/s.

    bwatflx = bmlt * dx * dy

    ! Compute total input of meltwater (m^3/s)
    total_flux_in = sum(bwatflx)  ! need global sum for parallel code
    if (verbose_bwat .and. main_task) then
       print*, ' '
       print*, 'total input basal melt flux (m^3/s):', total_flux_in
    endif

    flux_unrouted = 0.0d0

    ! Begin loop over points, highest first
    !TODO: need to parallelize this loop somehow

    do k = nn,1,-1

       ! Get x and y indices of current point
       i = sorted(k,1)
       j = sorted(k,2)

       ! If the flux to this cell is nonzero, then route it to adjacent downhill cells
       if (bwat_mask(i,j) == 1 .and. bwatflx(i,j) > 0.0d0) then

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

          !WHL - debug
          if (this_rank == rtest .and. i == itest .and. j == jtest) then
             print*, ' '
             print*, 'slope: task, i, j =', rtest, i, j
             print*, slope(:,1)
             print*, slope(:,0)
             print*, slope(:,-1)
             print*, 'sum(slope) =', sum(slope)
          endif

          ! If there are places for the water to drain, distribute it according to the flux-routing scheme:
          !  to the lowest-elevation neighbor (D8), the two lowest-elevation neighbors (Dinf), or
          !  all lower-elevation neighbors (FD8).
          ! The D8 and FD8 schemes have been tested with a simple dome problem.
          ! Dinf is less suited for the dome problem because there are many ties for 2nd greatest slope,
          !  so i2 and j2 for slope2 are not well defined.
          ! Note that the flux in the source cell is not zeroed.

          if (flux_routing_scheme == HO_FLUX_ROUTING_D8) then

             ! route to the adjacent cell with the lowest elevation
             slope1 = 0.0d0
             if (sum(slope) > 0.d0) then
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
                bwatflx(i1,j1) = bwatflx(i1,j1) + bwatflx(i,j)
             else
                flux_unrouted = flux_unrouted + bwatflx(i,j)
                print*, 'Warning: Cell with no downhill neighbors, i, j, bwatflx =', &
                     i, j, bwatflx(i,j)
             endif

             if (this_rank == rtest .and. i == itest .and. j == jtest) then
                print*, 'i1, j1, slope1 =', i1, j1, slope1
             endif

          !TODO - Remove Dinf scheme?
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
                      j2 = itmp
                   elseif (slope(ii,jj) > slope2) then
                      slope2 = slope(ii,jj)
                      i2 = ip
                      j2 = jp
                   endif
                enddo
             enddo

             sum_slope = slope1 + slope2
             if (sum_slope > 0.0d0) then    ! divide the flux between cells (i1,j1) and (i2,j2)
                if (slope1 > 0.0d0) then
                   bwatflx(i1,j1) = bwatflx(i1,j1) + bwatflx(i,j)*slope1/sum_slope
                endif
                if (slope2 > 0.0d0) then
                   bwatflx(i2,j2) = bwatflx(i2,j2) + bwatflx(i,j)*slope2/sum_slope
                endif
             else
                print*, 'Warning: Cell with no downhill neighbors, i, j =', i, j
             endif

             if (this_rank == rtest .and. i == itest .and. j == jtest) then
                print*, 'i1, j1, slope1:', i1, j1, slope1
                print*, 'i2, j2, slope2:', i2, j2, slope2
             endif

          elseif (flux_routing_scheme == HO_FLUX_ROUTING_FD8) then

             ! route to all adjacent downhill cells in proportion to grad(head)
             if (sum(slope) > 0.d0) then
                slope = slope / sum(slope)
                do jj = -1,1
                   do ii = -1,1
                      ip = i + ii
                      jp = j + jj
                      if (slope(ii,jj) > 0.d0) then
                         bwatflx(ip,jp) = bwatflx(ip,jp) + bwatflx(i,j)*slope(ii,jj)
                      endif
                   enddo
                enddo
             endif  ! sum(slope) > 0

          endif   ! flux_routing_scheme

       endif  ! bwat_mask = 1 and bwatflx > 0

    enddo  ! loop from high to low

    ! Identify cells just beyond the ice sheet margin, which can receive from upstream but not send downstream
    margin_mask = 0
    do j = 1, ny
       do i = 1, nx
          if (bwat_mask(i,j) == 0 .and. bwatflx(i,j) > 0.0d0) then
             margin_mask(i,j) = 1
          endif
       enddo
    enddo

    ! Compute total output of meltwater (m^3/s)

    !WHL - debug
!    print*, ' '
!    print*, 'Margin cells: i, j, bwatflx:'
    total_flux_out = 0.0d0
    do j = 1, ny
       do i = 1, nx
          if (margin_mask(i,j) == 1) then
             total_flux_out = total_flux_out + bwatflx(i,j)
          endif
       enddo
    enddo

    if (verbose_bwat .and. main_task) then
       print*, ' '
       print*, 'total output basal melt flux (m^3/s):', total_flux_out
       print*, 'total unrouted flux (m^3/s):', flux_unrouted
       print*, 'Sum:', total_flux_out + flux_unrouted
    endif

    !TODO - Add a bug check; should be equal

    ! clean up
    deallocate(sorted)

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

    use glimmer_physcon, only : grav
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
    p = 5
    if (verbose_bwat .and. this_rank == rtest) then
       print*, ' '
       print*, 'grad(head):'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(e10.3)',advance='no') grad_head(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

    p_exponent = 1.d0 / (p_flux_to_depth + 1.d0)

    ! Note: In Sommers et al. (2018), Eq. 6, the basal water flux q (m^2/s) is
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
       nx,   ny,    &
       phi,         &
       phi_mask)

    ! Fill depressions in the input field phi

    implicit none

    ! Input/output variables

    integer, intent(in) ::  &
         nx, ny              ! number of grid cells in each direction

    real(dp), dimension(nx,ny), intent(inout) :: &
         phi                 ! input field with depressions to be filled

    integer, dimension(nx,ny), intent(in) ::  &
         phi_mask            ! = 1 in the domain where depressions need to be filled, else = 0
                             ! corresponds to the grounded ice sheet for the flux-routing problem

    ! Local variables --------------------------------------

    real(dp), dimension(nx,ny) ::  &
         old_phi,          & ! old value of phi
         pool                ! identifies cells that need to be filled

    real(dp) :: pvs(9), max_val

    real(dp), parameter :: null = 1.d+20   ! large number
    integer :: flag, i, j

    integer :: count
    integer, parameter :: count_max = 200

!!    logical, parameter :: verbose_depressions = .false.
    logical, parameter :: verbose_depressions = .true.


    ! initialize

    flag = 1
    count = 0

    do while (flag == 1)

       count = count + 1
       if (verbose_depressions .and. main_task) then
          print*, ' '
          print*, 'fill_depressions, count =', count
       endif

       flag = 0
       old_phi = phi

       do j = 2, ny-1
          do i = 2, nx-1
             if (phi_mask(i,j) == 1) then

                if (any(old_phi(i-1:i+1,j-1:j+1) < old_phi(i,j))) then
                   pool(i,j) = 0
                else
                   pool(i,j) = 1
                end if

                if (pool(i,j) == 1) then
                   flag = 1
                   pvs = (/ old_phi(i-1:i+1,j-1), old_phi(i-1:i+1,j+1), old_phi(i-1:i+1,j) /)

                   where (pvs == old_phi(i,j))  ! equal to the original phi
                      pvs = null
                   end where

                   max_val = minval(pvs)

                   if (max_val /= null) then
                      phi(i,j) = max_val
                   else
                      flag = 0
                   end if

                   if (verbose_depressions) then
                      print*, 'flag, i, j, old phi, new phi:', flag, i, j, old_phi(i,j), phi(i,j)
                   endif

                end if   ! pool = 1

             end if   ! phi_mask = 1
          end do   ! i
       end do   ! j

       if (count > count_max) then
          call write_log('Hydrology error: too many iterations in fill_depressions', GM_FATAL)
       endif

    end do   ! flag = 1

  end subroutine fill_depressions

!==============================================================

  subroutine fix_flats(&
       nx,    ny,            &
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

    implicit none

    ! Input/output variables

    integer, intent(in) ::  &
         nx, ny,               & ! number of grid cells in each direction
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), dimension(nx,ny), intent(inout) :: &
         phi                    ! input field with flat regions to be fixed

    integer, dimension(nx,ny), intent(in) ::  &
         phi_mask               ! = 1 where any flat regions of phi will need to be fixed, else = 0
                                ! corresponds to the grounded ice sheet (bmlt_mask) for the flux-routing problem

    ! Local variables --------------------------------------

    real(dp), dimension(nx,ny) ::  &
         phi_input,           & ! input value of phi, before any corrections
         phi_new,             & ! new value of phi, after incremental corrections
         dphi1,               & ! sum of increments applied in step 1
         dphi2                  ! sum of increments applied in step 2

    integer, dimension(nx,ny) :: &
         flat_mask,           & ! = 1 for cells with phi_mask = 1 and without a downslope gradient, else = 0
         flat_mask_input,     & ! flat_mask as computed from phi_input
         n_uphill,            & ! number of uphill neighbors for each cell, as computed from input phi
         n_downhill             ! number of downhill neighbors for each cell, as computed from input phi

    logical, dimension(nx,ny) :: &
         incremented,         & ! = T for cells that have already been incremented (in step 2)
         incremented_neighbor   ! = T for cells that have not been incremented, but have an incremented neighbor
 
    logical :: finished         ! true when a loop has finished

    real(dp), parameter :: &
         phi_increment = 2.0d-5    ! fractional increment in phi (Garbrecht & Martz use 2.0e-5)

    integer :: i, j, ii, jj, ip, jp, p
    integer :: count
    integer, parameter :: count_max = 50

    !WHL - debug
!!    logical, parameter :: verbose_fix_flats = .false.
    logical, parameter :: verbose_fix_flats = .true.

    p = pdiag

    if (verbose_fix_flats .and. this_rank == rtest) then
       print*, ' '
       print*, 'In fix_flats, rtest, itest, jtest =', rtest, itest, jtest
       print*, ' '
       print*, 'input phi:'
       write(6,'(a3)',advance='no') '   '
       do i = itest-p, itest+p
          write(6,'(i10)',advance='no') i
       enddo
       write(6,*) ' '
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.5)',advance='no') phi(i,j)
          enddo
          write(6,*) ' '
       enddo
       write(6,*) ' '
       print*, 'mask:'
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(i10)',advance='no') phi_mask(i,j)
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

    call find_flats(&
         nx,    ny,           &
         itest, jtest, rtest, &
         phi_input,           &
         phi_mask,            &
         flat_mask_input)

    if (verbose_fix_flats .and. this_rank == rtest) then
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
    ! Repeat until all cells have a downslope gradient.

    do while(.not.finished)

       count = count + 1
       if (verbose_fix_flats .and. this_rank == rtest) then
          print*, ' '
          print*, 'step 1, count =', count
       endif

       do j = 2, ny-1
          do i = 2, nx-1
             if (flat_mask(i,j) == 1) then
                dphi1(i,j) = dphi1(i,j) + phi_increment
             endif
          enddo
       enddo

       if (verbose_fix_flats .and. this_rank == rtest) then
          print*, ' '
          print*, 'Updated dphi1/phi_increment:'
          do j = jtest+p, jtest-p, -1
             write(6,'(i6)',advance='no') j
             do i = itest-p, itest+p
                write(6,'(f10.1)',advance='no') dphi1(i,j)/ phi_increment
             enddo
             write(6,*) ' '
          enddo
       endif

       ! From the original flat region, identify cells that still have no downslope gradient.

       phi_new = phi_input + dphi1

       call find_flats(&
            nx,    ny,             &
            itest, jtest, rtest,   &
            phi_new,               &
            flat_mask_input,       &
            flat_mask)

!       call parallel_halo(flat_mask, parallel)

       if (verbose_fix_flats .and. this_rank == rtest) then
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

       ! If any flat cells remain, then repeat; else exit
       finished = .true.
       if (sum(flat_mask) > 0) then
          finished = .false.
       endif

       if (count > count_max) then
          call write_log('Hydrology error: abort in step 1 of fix_flats', GM_FATAL)
       endif

    enddo   ! step 1 finished

    ! Step 2: Gradient away from higher terrain

    dphi2 = 0.0d0
    incremented = .false.
    finished = .false.
    count = 0

    ! In the first pass, increment the elevation in all cells of the input flat region that are
    !  adjacent to higher terrain and have no adjacent downhill cell.

    do j = 2, ny-1
       do i = 2, nx-1
          if (flat_mask_input(i,j) == 1) then
             if (n_uphill(i,j) >= 1 .and. n_downhill(i,j) == 0) then
                dphi2(i,j) = dphi2(i,j) + phi_increment
                incremented(i,j) = .true.
             endif
          endif
       enddo
    enddo

!    call parallel_halo(incremented, parallel)

    if (verbose_fix_flats .and. this_rank == rtest) then
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

    ! If no cells are incremented in the first pass, then skip step 2.
    ! This will be the case if the flat region lies at the highest elevation in the domain.

    if (.not.any(incremented)) then
       if (verbose_fix_flats .and. this_rank == rtest) then
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
       if (verbose_fix_flats .and. this_rank == rtest) then
          print*, ' '
          print*, 'step 2, count =', count
       endif

       ! Identify cells that have not been incremented, but are adjacent to incremented cells
       incremented_neighbor = .false.
       do j = 2, ny-1
          do i = 2, nx-1
             if (flat_mask_input(i,j) == 1 .and. .not.incremented(i,j)) then
                do jj = -1,1
                   do ii = -1,1
                      ! If this is the centre point, ignore
                      if (ii == 0 .and. jj == 0) then
                         continue
                      else  ! check for an incremented neighbor
                         ip = i + ii
                         jp = j + jj
                         if (incremented(ip,jp)) then
                            incremented_neighbor(i,j) = .true.
                         endif
                      endif
                   enddo   ! ii
                enddo   ! jj
             endif   ! flat_mask = 1 and incremented = F
          enddo   ! i
       enddo   ! j

!       call parallel_halo(incremended_neighbor, parallel)

       ! Increment cells of type (1) and (2)
       ! Note: n_downhill was computed before step 1.

       do j = 2, ny-1
          do i = 2, nx-1
             if (incremented(i,j)) then
                dphi2(i,j) = dphi2(i,j) + phi_increment
             elseif (n_downhill(i,j) == 0 .and. incremented_neighbor(i,j)) then
                dphi2(i,j) = dphi2(i,j) + phi_increment
                incremented(i,j) = .true.
             endif
          enddo
       enddo

!       call parallel_halo(incremented, parallel)

       if (verbose_fix_flats .and. this_rank == rtest) then
          print*, ' '
          print*, 'incremented_neighbor:'
          do j = jtest+p, jtest-p, -1
             write(6,'(i6)',advance='no') j
             do i = itest-p, itest+p
                write(6,'(L10)',advance='no') incremented_neighbor(i,j)
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

       ! Check for cells that are in the input flat region and have not been incremented.
       ! If there are no such cells, then exit the loop.
       finished = .true.
       do j = 2, ny-1
          do i = 2, nx-1
             if (flat_mask_input(i,j) == 1 .and. .not.incremented(i,j)) then
                finished = .false.
                exit
             endif
          enddo
       enddo

       if (count > count_max) then
          call write_log('Hydrology error: abort in step 2 of fix_flats', GM_FATAL)
       endif

    enddo   ! step 2 finished


    ! Step 3:

    ! Add the increments from steps 1 and 2
    ! The result is a surface with gradients both toward lower terrain and away from higher terrain.

    phi = phi + dphi1 + dphi2

    ! Check for cells with flat_mask = 1 (no downslope gradient).
    ! Such cells are possible because of cancelling dphi1 and dphi2.

    count = 0
    finished = .false.

    do while (.not.finished)

       count = count + 1
       if (verbose_fix_flats .and. this_rank == rtest) then
          print*, ' '
          print*, 'step 3, count =', count
       endif

       ! Identify cells without downslope neighbors

       call find_flats(&
            nx,    ny,           &
            itest, jtest, rtest, &
            phi,                 &
            phi_mask,            &
            flat_mask)

       ! Add a half increment to any cells without downslope neighbors.
       ! If all cells have downslope neighbors, then exit.

       if (verbose_fix_flats .and. this_rank == rtest) then
          print*, 'sum(flat_mask) =', sum(flat_mask)
       endif

       if (sum(flat_mask) > 0) then
          where (flat_mask == 1)
             phi = phi + 0.5d0 * phi_increment
          endwhere
          finished = .false.
       else
          finished = .true.
       endif

       if (count > count_max) then
          call write_log('Hydrology error: abort in step 3 of fix_flats', GM_FATAL)
       endif

    enddo   ! step 3 finished

    if (verbose_fix_flats .and. this_rank == rtest) then
       print*, ' '
       print*, 'Final phi:'
       do j = jtest+p, jtest-p, -1
          write(6,'(i6)',advance='no') j
          do i = itest-p, itest+p
             write(6,'(f10.5)',advance='no') phi(i,j)
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
    ! Note: This definition includes some cells that have the same elevation as
    !       adjacent cells in the flat region, but have a nonzero downslope gradient.

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

  subroutine heights_sort(&
       nx,    ny,            &
       itest, jtest, rtest,  &
       head,  sorted)

    ! Create an array with the x and y location of each cell, sorted from from low to high values of head.
    ! TODO: Adapt for parallel code.  Sort only the locally owned grid cells?

    ! Input/output arguments

    integer, intent(in) ::  &
         nx, ny,               & ! number of grid cells in each direction
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), dimension(nx,ny), intent(in) :: &
         head                    ! hydraulic head (m), to be sorted from low to high

    integer, dimension(nx*ny,2), intent(inout) :: &
         sorted                  ! i and j indices of each cell, sorted from from low to high head

    ! Local variables

    integer :: nn, i, j, k
    real(dp), dimension(nx*ny) :: vect
    integer, dimension(nx*ny) :: ind

    nn = nx*ny

    ! Fill a work vector with head values
    k = 1
    do i = 1, nx
       do j = 1, ny
          vect(k) = head(i,j)
          k = k + 1
       enddo
    enddo

    ! Sort the vector from low to high values
    call indexx(vect, ind)

    ! Fill the 'sorted' array with the i and j values of each cell
    do k = 1, nn
       sorted(k,1) = floor(real(ind(k)-1)/real(ny)) + 1
       sorted(k,2) = mod(ind(k)-1,ny)+1
    enddo

    ! Fill the 'vect' array with head values
    ! Note: This array is not an output field; used only for a bug check

    do k = 1, nn
       vect(k) = head(sorted(k,1), sorted(k,2))
    enddo

    !WHL - debug
    if (verbose_bwat .and. this_rank == rtest) then
!!       print*, ' '
!!       print*, 'k, x, y, head:'
       do k = nn-20, nn
          vect(k) = head(sorted(k,1), sorted(k,2))
!!          print*, k, sorted(k,1), sorted(k,2), vect(k)
       enddo
    endif

  end subroutine heights_sort

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

    use glimmer_log

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
