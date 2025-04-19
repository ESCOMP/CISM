!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_transport.F90 - part of the Community Ice Sheet Model (CISM)  
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
!
! This module contains drivers for incremental remapping and upwind ice transport.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This version was created from ice_transport_driver in CICE, revision 313, 6 Jan. 2011.
! The repository is here: http://oceans11.lanl.gov/svn/CICE

  module glissade_transport

    use glimmer_global, only: dp
    use glimmer_log
    use glissade_remap, only: glissade_horizontal_remap, make_remap_mask, puny
    use cism_parallel, only: this_rank, main_task, nhalo, lhalo, uhalo, staggered_lhalo, staggered_uhalo, &
         parallel_type, parallel_reduce_max, parallel_reduce_sum, parallel_reduce_minloc, &
         parallel_globalindex, broadcast

    implicit none
    save
    private

    public :: glissade_transport_driver, glissade_check_cfl, &
              glissade_transport_setup_tracers, glissade_transport_finish_tracers

    public :: glissade_global_conservation, glissade_sum_mass_and_tracers, glissade_vertical_remap

    logical, parameter ::  &
         prescribed_area = .false.  ! if true, prescribe the area fluxed across each edge

    logical, parameter ::     &
         conservation_check = .true. ! if true, check global conservation

    logical, parameter :: verbose_ice_age = .false.

!=======================================================================

  contains

!=======================================================================
!
    subroutine glissade_transport_setup_tracers(model,        &
                                                transport_tracers_in)
      
      ! This subroutine copies all the 3D tracer fields into a single array for transport.

      use glide_types

      type(glide_global_type), intent(inout) :: model   ! derived type holding ice-sheet info

      logical, intent(in), optional ::  &
           transport_tracers_in

      integer :: nlyr    ! number of layers = upn-1

      integer :: k, nt

      logical :: transport_tracers

      integer :: i, j, nx, ny

      nx = model%general%ewn
      ny = model%general%nsn
      nlyr = model%general%upn - 1

      if (present(transport_tracers_in)) then
         transport_tracers = transport_tracers_in
      else
         transport_tracers = .true.  ! transport tracers by default
      endif

      if (.not.associated(model%geometry%tracers)) then

         ! first call for this instance; need to count tracers and allocate the tracer array

         model%geometry%ntracers = 0

         if (model%options%whichtemp == TEMP_PROGNOSTIC .or.  &
             model%options%whichtemp == TEMP_ENTHALPY) then
            model%geometry%ntracers = model%geometry%ntracers + 1
         endif

         if (model%options%which_ho_ice_age == HO_ICE_AGE_COMPUTE) then
            model%geometry%ntracers = model%geometry%ntracers + 1
         endif

         if (model%options%whichcalving == CALVING_DAMAGE) then
            model%geometry%ntracers = model%geometry%ntracers + 1
         endif

         ! add more tracers here, if desired

         ! make sure there is at least one tracer
         model%geometry%ntracers = max(model%geometry%ntracers,1)

         ! allocate and initialize tracer arrays
         allocate(model%geometry%tracers(nx,ny,model%geometry%ntracers,nlyr))
         allocate(model%geometry%tracers_usrf(nx,ny,model%geometry%ntracers))
         allocate(model%geometry%tracers_lsrf(nx,ny,model%geometry%ntracers))         

         model%geometry%tracers(:,:,:,:) = 0.0d0
         model%geometry%tracers_usrf(:,:,:) = 0.0d0
         model%geometry%tracers_lsrf(:,:,:) = 0.0d0

      endif   ! not(associated(model%geometry%tracers))

      ! fill the tracer array

      if (transport_tracers) then

         nt = 0

         ! start with temperature/enthalpy
         ! Note: temp/enthalpy values at upper surface (k=0) and lower surface (k=upn) are not transported,
         !       but these values are applied to new accumulation at either surface
         !       (in subroutine add_surface_and_basal_mass_balance)
         ! TODO: Set tracers_usrf to min(artm, 0.0) instead of temp(0) for the case that a cell
         !       is currently ice-free with temp(0) = 0?

         if (model%options%whichtemp == TEMP_PROGNOSTIC) then

            nt = nt + 1
            do k = 1, nlyr
               model%geometry%tracers(:,:,nt,k) = model%temper%temp(k,:,:)
            enddo
            model%geometry%tracers_usrf(:,:,nt) = model%temper%temp(0,:,:)
            model%geometry%tracers_lsrf(:,:,nt) = model%temper%temp(nlyr+1,:,:)

         elseif (model%options%whichtemp == TEMP_ENTHALPY) then

            nt = nt + 1
            do k = 1, nlyr
               model%geometry%tracers(:,:,nt,k) = model%temper%enthalpy(k,:,:)
            enddo
            model%geometry%tracers_usrf(:,:,nt) = model%temper%enthalpy(0,:,:)
            model%geometry%tracers_lsrf(:,:,nt) = model%temper%enthalpy(nlyr+1,:,:)

         endif

         ! damage parameter for prognostic calving scheme
         if (model%options%whichcalving == CALVING_DAMAGE) then

            nt = nt + 1
            do k = 1, nlyr
               model%geometry%tracers(:,:,nt,k) = model%calving%damage(k,:,:)
            enddo

            !Note: The damage for surface accumulation is set to the damage in layer 1,
            !       and the damage for basal freeze-on is set to the damage in layer nlyr.
            !      In other words, if the ice is damaged, then new snow/ice will not heal it,
            !       and if the ice is undamaged, new snow/ice will not increase the damage.
            !      This approach was suggested by Jeremy Bassis (8/4/15). Other hypotheses may be tested later.

            model%geometry%tracers_usrf(:,:,nt) = model%geometry%tracers(:,:,nt,1)
            model%geometry%tracers_lsrf(:,:,nt) = model%geometry%tracers(:,:,nt,nlyr)

         endif
         
         ! ice age parameter
         if (model%options%which_ho_ice_age == HO_ICE_AGE_COMPUTE) then

            nt = nt + 1
            do k = 1, nlyr
               model%geometry%tracers(:,:,nt,k) = model%geometry%ice_age(k,:,:)
            enddo

            ! Note: Upper surface accumulation has age = 0.
            !       Basal freeze-on, however, is assigned the same age as the bottom layer.
            !        This is based partly on the fact that it may originate from old melted ice.
            !        Also, the 2nd order vertical remapping constructs a gradient based on this value.
            !        A lower surface value of zero will likely lead to a poor approximation of the gradient.
            
            model%geometry%tracers_usrf(:,:,nt) = 0.0d0
            model%geometry%tracers_lsrf(:,:,nt) = model%geometry%tracers(:,:,nt,nlyr)
            
         endif

         ! add more tracers here, if desired

      endif  ! transport_tracers

    end subroutine glissade_transport_setup_tracers

!=======================================================================

    subroutine glissade_transport_finish_tracers(model)

      ! This subroutine copies the 3D tracer fields from a single tracer transport array
      !  back to individual tracer arrays.
      ! It also does halo updates for tracers.

      use glide_types

      type(glide_global_type), intent(inout) :: model   ! derived type holding ice-sheet info

      integer :: k, nt
      integer :: nlyr ! number of vertical layers = upn-1

      !WHL - debug
      integer :: i, j, nx, ny
      nx = model%general%ewn
      ny = model%general%nsn

      nlyr = model%general%upn - 1

      nt = 0

      ! start with temperature or enthapy, depending on model%options%whichtemp

      if (model%options%whichtemp == TEMP_PROGNOSTIC) then
         nt = nt+1
         do k = 1, nlyr
            model%temper%temp(k,:,:) = model%geometry%tracers(:,:,nt,k)
         enddo
      elseif (model%options%whichtemp == TEMP_ENTHALPY) then
         nt = nt+1
         do k = 1, nlyr
            model%temper%enthalpy(k,:,:) = model%geometry%tracers(:,:,nt,k)
         enddo
      endif

      ! damage parameter for prognostic calving scheme
      if (model%options%whichcalving == CALVING_DAMAGE) then
         nt = nt + 1
         do k = 1, nlyr
            model%calving%damage(k,:,:) = model%geometry%tracers(:,:,nt,k) 
         enddo
      endif

      ! ice age parameter
      if (model%options%which_ho_ice_age == HO_ICE_AGE_COMPUTE) then
         nt = nt + 1
         do k = 1, nlyr
            model%geometry%ice_age(k,:,:) = model%geometry%tracers(:,:,nt,k)
         enddo
      endif
      
      ! add more tracers here, if desired

    end subroutine glissade_transport_finish_tracers

!=======================================================================

    subroutine glissade_transport_driver(dt,                         &
                                         dx,           dy,           &
                                         nx,           ny,           &
                                         nlyr,         sigma,        &
                                         parallel,                   &
                                         itest, jtest, rtest,        &
                                         uvel,         vvel,         &
                                         thck,                       &
                                         ntracers,     tracers,      &
                                         tracers_usrf, tracers_lsrf, &
                                         vert_remap_accuracy,        &
                                         upwind_transport_in)

      ! This subroutine solves the transport equations for one timestep
      ! using the conservative remapping scheme developed by John Dukowicz
      ! and John Baumgardner and modified for sea ice by William
      ! Lipscomb and Elizabeth Hunke.
      !
      ! This scheme preserves monotonicity of ice area and tracers.  That is,
      ! it does not produce new extrema.  It is second-order accurate in space,
      ! except where gradients are limited to preserve monotonicity. 
      !
      ! Optionally, the remapping scheme can be replaced with a simple
      ! first-order upwind scheme.
      !
      ! author William H. Lipscomb, LANL
      !
      ! input/output arguments

      real(dp), intent(in) ::  &
         dt,                   &! time step (s)
         dx, dy                 ! gridcell dimensions (m)
                                ! (cells assumed to be rectangular)

      integer, intent(in) ::   &
         nx, ny,               &! horizontal array size
         nlyr                   ! number of vertical layers

      real(dp), dimension(nlyr+1), intent(in) ::  &
         sigma                  ! layer interfaces in sigma coordinates
                                ! top sfc = 0, bottom sfc = 1

      type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

      integer, intent(in) ::  &
         itest, jtest, rtest    ! coordinates of diagnostic point

      real(dp), dimension(nlyr+1,nx-1,ny-1), intent(in) ::  &
         uvel, vvel             ! horizontal velocity components (m/s)
                                ! (defined at horiz cell corners, vertical interfaces)

      real(dp), dimension(nx,ny), intent(inout) ::  &
         thck                   ! ice thickness (m), defined at horiz cell centers

      integer, intent(in) ::  &
         ntracers               ! number of tracers to be transported

      !TODO - Make the tracer arrays optional arguments?
      real(dp), dimension(nx,ny,ntracers,nlyr), intent(inout) ::  &
         tracers                ! set of 3D tracer arrays, packed into a 4D array

      real(dp), dimension(nx,ny,ntracers), intent(in) :: &
         tracers_usrf,         &! tracer values associated with accumulation at upper surface
         tracers_lsrf           ! tracer values associated with freeze-on at lower surface

      integer, intent(in) ::  &
         vert_remap_accuracy    ! order of accuracy for vertical remapping
                                ! HO_VERTICAL_REMAP_FIRST_ORDER or HO_VERTICAL_REMAP_SECOMD_ORDER

      logical, intent(in), optional ::  &
         upwind_transport_in    ! if true, do first-order upwind transport

      ! local variables

      integer ::     &
         i, j, k         ,&! cell indices
         ilo,ihi,jlo,jhi ,&! beginning and end of physical domain
         nt                ! tracer index

      real(dp), dimension (nx,ny) ::     &
         thck_mask         ! = 1. if ice is present, = 0. otherwise

      real(dp), dimension (nx-1,ny-1) ::     &
         uvel_layer      ,&! uvel averaged to layer midpoint (m/s)
         vvel_layer        ! vvel averaged to layer midpoint (m/s)

      real(dp), dimension (nx,ny,nlyr) ::     &
         thck_layer        ! ice layer thickness (m)

      integer ::     &
         icells            ! number of cells with ice

      integer, dimension(nx*ny) ::     &
         indxi, indxj      ! compressed i/j indices

      !-------------------------------------------------------------------
      ! If prescribed_area is true, the area of each departure region is
      !  computed in advance (e.g., by taking the divergence of the
      !  velocity field and passed to locate_triangles.  The departure
      !  regions are adjusted to obtain the desired area.
      ! If false, edgearea is computed in locate_triangles and passed out.
      !-------------------------------------------------------------------

      real(dp), dimension(nx,ny) ::   &
         edgearea_e     ,&! area of departure regions for east edges
         edgearea_n       ! area of departure regions for north edges

      real(dp) ::     &
         msum_init,      &! initial global ice mass
         msum_final       ! final global ice mass

      real(dp), dimension(ntracers) ::     &
         mtsum_init,     &! initial global ice mass*tracer
         mtsum_final      ! final global ice mass*tracer

      logical ::     &
         errflag          ! true if energy is not conserved

      character(len=100) :: message

      real(dp), dimension (:,:,:), allocatable :: &
         worku            ! work array

      real(dp), dimension(nx,ny) ::      &
         uee, vnn            ! cell edge velocities for upwind transport

      logical ::     &
         upwind_transport    ! if true, do first-order upwind transport

      !-------------------------------------------------------------------
      ! Initialize
      !-------------------------------------------------------------------

      if (present(upwind_transport_in)) then
         upwind_transport = upwind_transport_in
      else
         upwind_transport = .false.
      endif

      errflag = .false.

      !Note: (ilo,ihi) and (jlo,jhi) are the lower and upper bounds of the local domain
      ! (i.e., grid cells owned by this processor).

      ilo = nhalo + 1
      ihi = nx - nhalo
      jlo = nhalo + 1
      jhi = ny - nhalo

      !-------------------------------------------------------------------
      ! NOTE: Mass and tracer arrays (thck, temp, etc.) must be updated in 
      !       halo cells before this subroutine is called. 
      !-------------------------------------------------------------------

      !-------------------------------------------------------------------
      ! Fill layer thickness array.
      !-------------------------------------------------------------------

      do k = 1, nlyr
         thck_layer(:,:,k) = thck(:,:) * (sigma(k+1) - sigma(k))
      enddo

      !-------------------------------------------------------------------
      ! Compute initial values of globally conserved quantities (optional)
      !-------------------------------------------------------------------

      if (conservation_check) then

         call glissade_sum_mass_and_tracers(&
              nx,                ny,              &
              nlyr,              ntracers,        &
              thck_layer(:,:,:), msum_init,       &
              tracers(:,:,:,:),  mtsum_init(:))
      endif

      !-------------------------------------------------------------------
      ! Horizontal transport of ice thickness and tracers
      ! Two options:
      ! (1) First-order accurate upwind scheme
      ! (2) Second-order accurate incremental remapping scheme
      !-------------------------------------------------------------------

      if (upwind_transport) then

         allocate (worku(nx,ny,0:ntracers))

         do k = 1, nlyr

            ! Average corner velocities at layer interfaces to 
            ! edge velocities at layer midpoints.
      
            do j = jlo, jhi
            do i = ilo-1, ihi   ! include west edge of local domain
               uee(i,j) = 0.25d0 * (uvel(k,  i,j) + uvel(k,  i,j-1)    &
                                  + uvel(k+1,i,j) + uvel(k+1,i,j-1))
            enddo
            enddo


            do j = jlo-1, jhi   ! include south edge of local domain
            do i = ilo, ihi
               vnn(i,j) = 0.25d0 * (vvel(k,  i,j) + vvel(k,  i-1,j)    &
                                  + vvel(k+1,i,j) + vvel(k+1,i-1,j))
            enddo
            enddo

            ! Fill work array for transport

            worku(:,:,0) = thck_layer(:,:,k)
            do nt = 1, ntracers
               worku(:,:,nt) = thck_layer(:,:,k) * tracers(:,:,nt,k)
            enddo

            !-----------------------------------------------------------------
            ! Upwind transport
            !-----------------------------------------------------------------
 
            do nt = 0, ntracers
               call upwind_field (nx,             ny,                  &
                                  ilo, ihi,       jlo, jhi,            &
                                  dx,             dy,                  &
                                  dt,             worku(:,:,nt),       &
                                  uee(:,:),       vnn    (:,:))
            enddo   ! ntracers

            ! Recompute tracers

            thck_layer(:,:,k) = worku(:,:,0)
            do nt = 1, ntracers
               do j = jlo, jhi
               do i = ilo, ihi
                  if (thck_layer(i,j,k) > puny) then
                     tracers(i,j,nt,k) = worku(i,j,nt) / thck_layer(i,j,k)
                  else
                     tracers(i,j,nt,k) = 0.d0
                  endif
               enddo   ! i
               enddo   ! j
            enddo      ! ntracers

         enddo         ! nlyr

         deallocate (worku)

      else    ! remapping transport

      !-------------------------------------------------------------------
      ! Define a mask: = 1 where ice is present (thck > 0), = 0 otherwise         
      ! The mask is used to prevent tracer values in cells without ice from
      !  being used to compute tracer gradients.
      !-------------------------------------------------------------------

         call make_remap_mask (nx,           ny,                 &
                               ilo, ihi,     jlo, jhi,           &
                               nhalo,        icells,             &
                               indxi(:),     indxj(:),           &
                               thck(:,:),    thck_mask(:,:))

      !WHL - debug
!      k = 2
!      write(6,*) 'Before remapping, tracer, k =', k
!      do j = ny, 1, -1
!         write(6,'(i6)',advance='no') j
!         do i = 5, nx-5
!            write(6,'(f8.3)',advance='no') tracers(i,j,1,k)
!         enddo
!         write(6,*) ' '
!      enddo

      !-------------------------------------------------------------------
      ! Remap ice thickness and tracers; loop over layers
      !-------------------------------------------------------------------

         do k = 1, nlyr

            ! Average velocities to the midpoint of the layer

            uvel_layer(:,:) = 0.5d0 * (uvel(k,:,:) + uvel(k+1,:,:))
            vvel_layer(:,:) = 0.5d0 * (vvel(k,:,:) + vvel(k+1,:,:))

            edgearea_e(:,:) = 0.d0
            edgearea_n(:,:) = 0.d0

            ! If prescribed_area is true, compute edgearea by taking the divergence
            !  of the velocity field.

            if (prescribed_area) then

               do j = jlo, jhi
               do i = ilo-1, ihi
                  edgearea_e(i,j) = (uvel_layer(i,j) + uvel_layer(i,j-1))     &
                                    * 0.5d0 * dy * dt
               enddo
               enddo
 
               do j = jlo-1, jhi
               do i = ilo, ihi
                  edgearea_n(i,j) = (vvel_layer(i,j) + vvel_layer(i-1,j))     &
                                   * 0.5d0 * dx * dt
               enddo
               enddo

            endif

         !-------------------------------------------------------------------
         ! Main remapping routine: Step ice thickness and tracers forward in time.
         !-------------------------------------------------------------------

            call glissade_horizontal_remap (dt,                                  &
                                            dx,                dy,               &
                                            nx,                ny,               &
                                            ntracers,          nhalo,            &
                                            parallel,                            &
                                            thck_mask(:,:),    icells,           &
                                            indxi(:),          indxj(:),         &
                                            uvel_layer(:,:),   vvel_layer(:,:),  &
                                            thck_layer(:,:,k), tracers(:,:,:,k), &
                                            edgearea_e(:,:),   edgearea_n(:,:))

         enddo    ! nlyr

      endif       ! remapping v. upwind transport

      !-------------------------------------------------------------------
      ! Check that mass and mass*tracers are exactly conserved by transport.
      ! Note: Conservation errors will occur if the global domain is open
      !       and ice has left the domain. So depending on the application,
      !       there may or may not be a problem when ice is not conserved.
      !-------------------------------------------------------------------

      if (conservation_check) then

         ! Compute new values of globally conserved quantities.
         ! Assume gridcells of equal area, ice of uniform density.

         call glissade_sum_mass_and_tracers(&
              nx,                ny,              &
              nlyr,              ntracers,        &
              thck_layer(:,:,:), msum_final,      &
              tracers(:,:,:,:),  mtsum_final(:))

         ! Check conservation

         if (main_task) then

            call glissade_global_conservation(&
                 msum_init,     msum_final,      &
                 errflag,       0.0d0,           & ! melt_potential = 0 for this check
                 ntracers,                       &
                 mtsum_init,    mtsum_final)

            if (errflag) then
               write(message,*) 'WARNING: Conservation error in glissade_horizontal_remap'
!               call write_log(message,GM_FATAL)      ! uncomment to make conservation errors fatal
               call write_log(message,GM_DIAGNOSTIC)  ! uncomment for debugging
               write(message,*) 'May be OK if global domain is open'
               call write_log(message,GM_DIAGNOSTIC)  ! uncomment for debugging
            endif

         endif   ! main_task

      endif      ! conservation_check

      !-------------------------------------------------------------------
      ! Interpolate tracers back to sigma coordinates
      !-------------------------------------------------------------------

      call glissade_vertical_remap(&
           nx,             ny,       &
           nlyr,           ntracers, &
           parallel,                 &
           itest,  jtest,  rtest,    &
           sigma(:),                 &
           thck_layer(:,:,:),        &
           tracers(:,:,:,:),         &
           tracers_usrf(:,:,:),      &
           tracers_lsrf(:,:,:),      &
           vert_remap_accuracy)

      !-------------------------------------------------------------------
      ! Final conservation check: Check that mass and mass*tracers are
      !                           conserved by vertical remapping.
      !-------------------------------------------------------------------

      if (conservation_check) then

         ! Update msum_init and mtsum_init
         msum_init = msum_final           ! msum_final computed above after horizontal transport
         mtsum_init(:) = mtsum_final(:)   ! mtsum_final computed above after horizontal transport
 
         ! Compute new values of globally conserved quantities.
         ! Assume gridcells of equal area, ice of uniform density.

         call glissade_sum_mass_and_tracers(&
              nx,                ny,              &
              nlyr,              ntracers,        &
              thck_layer(:,:,:), msum_final,      &
              tracers(:,:,:,:),  mtsum_final(:))

         ! Check conservation

         if (main_task) then

            call glissade_global_conservation(&
                 msum_init,     msum_final,  &
                 errflag,       0.0d0,       &  ! melt potential = 0 for this check
                 ntracers,                   &
                 mtsum_init,    mtsum_final)

            if (errflag) then
               write(message,*) 'WARNING: Conservation error in glissade_vertical_remap'
!               call write_log(message,GM_FATAL)      ! uncomment to make conservation errors fatal
               call write_log(message,GM_DIAGNOSTIC)  ! uncomment for debugging
            endif

         endif   ! main_task

      endif      ! conservation_check

      !-------------------------------------------------------------------
      ! Recompute thickness
      ! Note: Halo updates for thickness and tracers are done in glissade_transport_solve.
      !-------------------------------------------------------------------

      thck(:,:) = 0.d0
      do k = 1, nlyr
         thck(:,:) = thck(:,:) + thck_layer(:,:,k)
      enddo

    end subroutine glissade_transport_driver

!=======================================================================

    subroutine glissade_check_cfl(ewn,     nsn,      nlyr,          & 
                                  dew,     dns,      sigma,         &
                                  parallel,                         &
                                  stagthk, dusrfdew, dusrfdns,      &
                                  uvel,    vvel,     deltat,        &
                                  adaptive_cfl_threshold,           &
                                  allowable_dt_adv,  allowable_dt_diff)

      ! Calculate maximum allowable time step based on both 
      ! advective and diffusive CFL limits.
      !
      ! author Matt Hoffman, LANL, March 2014
      !
      ! input/output arguments

      integer, intent(in) ::     &
         ewn, nsn    ! number of cells in the x, y dimensions

      integer, intent(in) ::     &
         nlyr        ! number of vertical layers (layer centers)

      real(dp), intent(in) :: &
         dew, dns    ! grid spacing in x, y (not assumed to be equal here), dimensional m

      real(dp), dimension(:), intent(in) :: &
         sigma       ! vertical coordinate spacing

      type(parallel_type), intent(in) :: &
         parallel    ! info for parallel communication

      real(dp), dimension(:,:), intent(in) :: &
         stagthk     ! thickness on the staggered grid, dimensional m

      real(dp), dimension(:,:), intent(in) :: &
         dusrfdew, dusrfdns    ! slope in x,y directions on the staggered grid, dimensionless m/m

      real(dp), dimension(:,:,:), intent(in) :: &
         uvel, vvel  ! 3-d x,y velocity components on the staggered grid, dimensional m/yr

      real(dp), intent(in) :: &
         deltat      ! model deltat (yrs)

      real(dp), intent(in) :: &
         adaptive_cfl_threshold  ! threshold for adaptive subcycling
                                 ! if = 0, there is no adaptive subcycling; code aborts when CFL > 1

      real(dp), intent(out) :: &
         allowable_dt_adv     ! maximum allowable dt (yrs) based on advective CFL 

      real(dp), intent(out) :: &
         allowable_dt_diff    ! maximum allowable dt (yrs) based on diffusive CFL 

      ! Local variables
      integer :: k
      integer :: xs, xe, ys, ye  ! start and end indices for locally owned cells on the staggered grid in the x and y directions
      real(dp), dimension(nlyr, ewn-1, nsn-1) :: uvel_layer, vvel_layer  ! velocities at layer midpoints, stag. grid
      real(dp), dimension(nlyr, ewn-1, nsn-1) :: flux_layer_ew, flux_layer_ns  ! flux for each layer, stag. grid
      real(dp), dimension(ewn-1, nsn-1) :: flux_ew, flux_ns  ! flux for entire thickness, stag. grid

      real(dp) :: maxuvel, maxvvel, maxvel ! maximum velocity in either direction and in both
      real(dp) :: allowable_dt_diff_here ! temporary calculation at each cell of allowable_dt_diff
      real(dp) :: my_allowable_dt_adv  ! allowable_dt_adv on this processor
      real(dp) :: my_allowable_dt_diff ! allowable_dt_diff on this processor
      integer :: i, j
      real(dp) :: slopemag  ! the magnitude of the surface slope
      real(dp) :: slopedirx, slopediry  ! the unit vector of the slope direction
      real(dp) :: flux_downslope  ! The component of the flux in the downslope direction
      integer :: ierr ! flag for CFL violation
      integer :: procnum  ! processor on which minimum allowable time step occurs
      integer, dimension(3) :: indices_adv  ! z,x,y indices (stag. grid) of where the min. allow. time step occurs for the advective CFL
      integer, dimension(2) :: indices_diff  ! x and y indices (stag. grid) of where the min. allow. time step occurs for  the  diffusive CFL
      integer, dimension(3) :: indices_adv_global  ! z,x,y indices (stag. grid) corresponding to indices_adv, but in global rather than local index space
      integer, dimension(2) :: indices_diff_global ! x,y indices (stag. grid) corresponding to indices_diff, but in global rather than local index space
      character(len=12)  :: dt_string, xpos_string, ypos_string
      character(len=300) :: message
      ierr = 0

      ! Setup some mesh information - start and end indices for locally owned cells on the staggered grid in the x and y directions
      xs = 1+staggered_lhalo
      xe = ewn - 1 - staggered_uhalo
      ys = 1+staggered_lhalo
      ye = nsn - 1 - staggered_uhalo

      ! ------------------------------------------------------------------------
      ! Advective CFL
      ! TODO use depth-averaged velocity or layer-by-layer (top layer only?), or something else (BB09)?
      ! For now check all layers

      ! Calculate depth-averaged flux and velocity on the B-grid.
      ! The IR code basically uses a B-grid, the FO-Upwind method uses a C-grid.
      ! The B-grid calculation should be more conservative because that is where the
      ! velocities are calculated so there will be no averaging.
      ! (Also, IR is the primary advection method, so make this check most appropriate for that.)

      do k = 1, nlyr
         ! Average velocities to the midpoint of the layer
         uvel_layer(k,:,:) = 0.5d0 * (uvel(k,:,:) + uvel(k+1,:,:))
         vvel_layer(k,:,:) = 0.5d0 * (vvel(k,:,:) + vvel(k+1,:,:))
         ! calculate flux components for this layer
         flux_layer_ew(k,:,:) = uvel_layer(k,:,:) * stagthk(:,:) * (sigma(k+1) - sigma(k))
         flux_layer_ns(k,:,:) = vvel_layer(k,:,:) * stagthk(:,:) * (sigma(k+1) - sigma(k))
      enddo
      flux_ew(:,:) = sum(flux_layer_ew, 1)
      flux_ns(:,:) = sum(flux_layer_ns, 1)

      ! Advective CFL calculation - using all layers. Check locally owned cells only!
      maxuvel = maxval(abs(uvel_layer(:,xs:xe,ys:ye)))
      maxvvel = maxval(abs(vvel_layer(:,xs:xe,ys:ye)))
      ! Determine in which direction the max velocity is - Assuming dx=dy here!
      if (maxuvel > maxvvel) then
!         print *, 'max vel is in uvel'
         maxvel = maxuvel
         indices_adv = maxloc(abs(uvel_layer(:,xs:xe,ys:ye)))
      else
!         print *, 'max vel is in vvel'
         maxvel = maxvvel
         indices_adv = maxloc(abs(vvel_layer(:,xs:xe,ys:ye)))
      endif
      indices_adv(2:3) = indices_adv(2:3) + staggered_lhalo  ! want the i,j coordinates WITH the halo present -
                                                             ! we got indices into the slice of owned cells
      ! Finally, determine maximum allowable time step based on advective CFL condition.
      my_allowable_dt_adv = dew / (maxvel + 1.0d-20)

      ! ------------------------------------------------------------------------
      ! Diffusive CFL
      ! Estimate diffusivity using the relation that the 2-d flux Q=-D grad h and Q=UH, 
      ! where h is surface elevation, D is diffusivity, U is 2-d velocity vector, and H is thickness
      ! Solving for D = UH/-grad h
      !TODO - Modify this loop to consider only grounded ice.  The diffusive CFL computed for floating ice
      !       usually is unnecessarily small for HO problems.

      my_allowable_dt_diff = 1.0d20  ! start with a huge value
      indices_diff(:) = 1 ! Initialize these to something, on the off-chance they never get set... (e.g., no ice on this processor)
      ! Loop over locally-owned cells only!
      do j = ys, ye
         do i = xs, xe
            if (stagthk(i,j) > 0.0d0) then  ! don't bother doing all this for non-ice cells
                ! Find downslope vector
                slopemag = dsqrt(dusrfdew(i,j)**2 + dusrfdns(i,j)**2 + 1.0d-20)
                slopedirx = dusrfdew(i,j) / slopemag
                slopediry = dusrfdns(i,j) / slopemag
                ! Estimate flux in the downslope direction (Flux /dot -slopedir)
                flux_downslope = flux_ew(i,j) * (-1.0d0) * slopedirx + flux_ns(i,j) * (-1.0d0) * slopediry  ! TODO check signs here - they seem ok
                !!! Estimate diffusivity in the downslope direction only
                !!diffu = flux_downslope / slopemag
                !!my_allowable_dt_diff = 0.5d0 * dew**2 / (diffu + 1.0e-20)  ! Note: assuming diffu is isotropic here.
                ! DCFL: dt = 0.5 * dx**2 / D = 0.5 * dx**2 * slopemag / flux_downslope
                allowable_dt_diff_here = 0.5d0 * dew**2 * slopemag / (flux_downslope + 1.0e-20)  ! Note: assuming diffu is isotropic here.  assuming dx=dy
                if (allowable_dt_diff_here < 0.0d0) allowable_dt_diff_here = 1.0d20 ! ignore negative dt's (upgradient flow due to membrane stresses)
                if (allowable_dt_diff_here < my_allowable_dt_diff) then
                   my_allowable_dt_diff = allowable_dt_diff_here
                   indices_diff(1) = i
                   indices_diff(2) = j
                endif
            endif
         enddo
      enddo

      ! Determine location limiting the DCFL
!      print *, 'diffu dt', my_allowable_dt_diff, indices_diff(1), indices_diff(2)

      ! Optional print of local limiting dt on each procesor
      !print *,'LOCAL ADV DT, POSITION', my_allowable_dt_adv, indices_adv(2), indices_adv(3)
      !print *,'LOCAL DIFF DT, POSITION', my_allowable_dt_diff, indices_diff(1), indices_diff(2)

      ! ------------------------------------------------------------------------
      ! Now check for errors

      ! Perform global reduce for advective time step and determine where in the domain it occurs
      call parallel_reduce_minloc(xin=my_allowable_dt_adv, xout=allowable_dt_adv, xprocout=procnum)

      if (deltat > allowable_dt_adv) then
          ierr = 1  ! advective CFL violation

          ! Get position of the limiting location - do this only if an error message is needed to avoid 2 MPI comms
          indices_adv_global(1) = indices_adv(1)
          call parallel_globalindex(indices_adv(2), indices_adv(3), &
                                    indices_adv_global(2), indices_adv_global(3),  &
                                    parallel)
          ! Note: This subroutine assumes the scalar grid, but should work fine for the stag grid too
          ! indices_adv_global now has i,j on the global grid for this proc's location
          call broadcast(indices_adv_global(2), proc=procnum)
          call broadcast(indices_adv_global(3), proc=procnum)
          ! indices_adv_global now has i,j on the global grid for the limiting proc's location

          write(dt_string,'(f12.6)') allowable_dt_adv
          write(xpos_string,'(i12)') indices_adv_global(2)
          write(ypos_string,'(i12)') indices_adv_global(3)
          write(message,*) 'Advective CFL violation!  Maximum allowable time step for advective CFL condition is ' &
               // trim(adjustl(dt_string)) // ' yr, limited by global position i=' &
               // trim(adjustl(xpos_string)) // ' j=' //trim(adjustl(ypos_string))
          call write_log(trim(message),GM_WARNING)

          ! If adaptive subcyling is allowed, then make the code abort for egregious CFL violations,
          ! (defined as deltat > 10 * allowable_dt_adv), to prevent excessive subcycling.

          if (main_task .and. adaptive_cfl_threshold > 0.0d0) then
             if (deltat > 10.d0 * allowable_dt_adv) then
                print*, 'deltat, allowable_dt_adv, ratio =', deltat, allowable_dt_adv, deltat/allowable_dt_adv
                call write_log('Aborting with CFL violation', GM_FATAL)
             endif
             !WHL - debug
             if (deltat > allowable_dt_adv) then
                print*, 'deltat, allowable_dt_adv, ratio =', deltat, allowable_dt_adv, deltat/allowable_dt_adv
                print*, '  Limited by position', indices_adv_global(2), indices_adv_global(3)
             endif
          endif

      endif

      ! Perform global reduce for diffusive time step and determine where in the domain it occurs
      call parallel_reduce_minloc(xin=my_allowable_dt_diff, xout=allowable_dt_diff, xprocout=procnum)

      if (deltat > allowable_dt_diff) then
          ! Get position of the limiting location - do this only if an error message is needed to avoid 2 MPI comms
          call parallel_globalindex(indices_diff(1), indices_diff(2), &
                                    indices_diff_global(1), indices_diff_global(2),  &
                                    parallel)
          ! Note: this subroutine assumes the scalar grid, but should work fine for the stag grid too
          ! indices_diff_global now has i,j on the global grid for this proc's location
          call broadcast(indices_diff_global(1), proc=procnum)
          call broadcast(indices_diff_global(2), proc=procnum)
          ! indices_diff_global now has i,j on the global grid for the limiting proc's location

          write(dt_string,'(f12.6)') allowable_dt_diff
          write(xpos_string,'(i12)') indices_diff_global(1)
          write(ypos_string,'(i12)') indices_diff_global(2)

          !WHL - Commenting out this warning for now, because the diffusive CFL violation is rarely meaningful for HO runs
!!          write(message,*) 'Diffusive CFL violation!  Maximum allowable time step for diffusive CFL condition is ' &
!!               // trim(adjustl(dt_string)) // ' yr, limited by global position i=' &
!!               // trim(adjustl(xpos_string)) // ' j=' //trim(adjustl(ypos_string))
!!          ! Diffusive CFL violation is just a warning (because it may be overly restrictive as currently formulated)
!!          call write_log(trim(message),GM_WARNING)    
!!          write(message,*) &
!!               '(Note the currently implemented diffusive CFL calculation may be overly restrictive for higher-order dycores.)'
!!          call write_log(trim(message))
      endif

      ! TODO enable this fatal error after more testing!
      ! Now that we have checked both, throw a fatal error for an ACFL violation
      !if (ierr == 1) then
      !   call write_log('Advective CFL violation is a fatal error.  See log for details.', GM_FATAL)
      !endif

    end subroutine glissade_check_cfl

!=======================================================================

    subroutine glissade_sum_mass_and_tracers(&
         nx,         ny,        &
         nlyr,       ntracer,   &
         thck_layer, msum,      &
         tracer,     mtsum)

      ! Compute values of globally conserved quantities.
      ! Assume gridcells of equal area (dx*dy), ice of uniform density.
      
      ! Input/output arguments

      integer, intent(in) ::   &
         nx, ny,               &! horizontal array size
         nlyr,                 &! number of vertical layers
         ntracer                ! number of tracers

      real(dp), dimension (nx,ny,nlyr), intent(in) ::     &
         thck_layer             ! ice layer thickness

      real(dp), intent(out) ::   &
         msum                   ! total mass (actually thickness, measured in m)

      real(dp), dimension (nx,ny,ntracer,nlyr), intent(in), optional ::   &
         tracer                 ! tracer values

      real(dp), dimension(ntracer), intent(out), optional ::   &
         mtsum                  ! total mass*tracer

      ! Local arguments

      integer :: i, j, nt

      msum  = 0.d0
      if (present(mtsum)) mtsum(:) = 0.d0

      do j = 1+nhalo, ny-nhalo
         do i = 1+nhalo, nx-nhalo
               
            ! accumulate ice mass and mass*tracers
            ! (actually, accumulate thickness, assuming rhoi*dx*dy is the same for each cell)

            msum = msum + sum(thck_layer(i,j,:))

            if (present(mtsum)) then
               do nt = 1, ntracer
                  mtsum(nt) =  mtsum(nt) + sum(tracer(i,j,nt,:)*thck_layer(i,j,:))
               enddo
            endif

         enddo    ! i
      enddo       ! j

      msum = parallel_reduce_sum(msum)
      if (present(mtsum)) mtsum = parallel_reduce_sum(mtsum)

    end subroutine glissade_sum_mass_and_tracers

!=======================================================================
!
    subroutine glissade_global_conservation(&
         msum_init,  msum_final,         &
         errflag,    melt_potential_in,  &
         ntracer,                        &
         mtsum_init, mtsum_final)

      ! Check whether values of conserved quantities have changed.
      ! An error probably means that ghost cells are treated incorrectly.
      !
      ! author William H. Lipscomb, LANL
      !
      ! input/output arguments

      real(dp), intent(in) ::     &
         msum_init   ,&! initial global ice mass 
         msum_final    ! final global ice mass

      logical, intent(out) ::     &
         errflag       ! true if there is a conservation error

      real(dp), intent(in), optional :: &
         melt_potential_in   ! total thickness (m) of additional ice that could be melted
                             ! by available acab/bmlt in columns that are completely melted

      integer, intent(in), optional ::  &
         ntracer       ! number of tracers

      real(dp), dimension(:), intent(in), optional :: &
         mtsum_init  ,&! initial global ice mass*tracer
         mtsum_final   ! final global ice mass*tracer

      character(len=100) :: message

      integer ::      &
         nt            ! tracer index

      real(dp) ::         &
         melt_potential,  &! melt_potential_in (if present), else = 0 
         diff              ! difference between initial and final values

      if (present(melt_potential_in)) then
         melt_potential = melt_potential_in
      else
         melt_potential = 0.d0
      endif

      errflag = .false.

      if (msum_init > puny) then
         diff = (msum_final - melt_potential) - msum_init
         if (abs(diff/msum_init) > puny) then
            errflag = .true.
            write (message,*) 'glissade_transport: ice mass conservation error'
            call write_log(message)
            write (message,*) 'Initial global mass =', msum_init
            call write_log(message)
            write (message,*) 'Final global mass =', msum_final
            call write_log(message)
            write (message,*) 'Melt potential =', melt_potential
            call write_log(message)
            write (message,*) 'Final global mass (adjusted for melt potential) =', msum_final - melt_potential
            call write_log(message)
            write (message,*) 'Absolute error =', diff
            call write_log(message)
            write (message,*) 'Fractional error =', abs(diff)/msum_init
            call write_log(message)
         endif
      endif

      if (present(mtsum_init)) then
         do nt = 1, ntracer
            if (abs(mtsum_init(nt)) > puny) then
               diff = mtsum_final(nt) - mtsum_init(nt)
               if (abs(diff/mtsum_init(nt)) > puny) then
                  errflag = .true.
                  write (message,*) 'glissade_transport: mass*tracer conservation error'
                  call write_log(message)
                  write (message,*) 'tracer index =', nt
                  call write_log(message)
                  write (message,*) 'Initial global mass*tracer =', mtsum_init(nt)
                  call write_log(message)
                  write (message,*) 'Final global mass*tracer =', mtsum_final(nt)
                  call write_log(message)
                  write (message,*) 'Fractional difference =', abs(diff)/mtsum_init(nt)
                  call write_log(message)
               endif
            endif
         enddo
      endif                     ! present(mtsum_init)

    end subroutine glissade_global_conservation

!----------------------------------------------------------------------

  subroutine glissade_vertical_remap(&
       nx,           ny,       &
       nlyr,         ntracer,  &
       parallel,               &
       itest, jtest, rtest,    &
       sigma,                  &
       hlyr,      trcr,        &
       trcr_usrf, trcr_lsrf,   &
       vert_remap_accuracy)

    use cism_parallel, only: parallel_reduce_max, parallel_reduce_min
    use glimmer_paramets, only: tim0
    use glimmer_physcon, only: scyr

    ! Conservative remapping of tracer fields from one set of vertical 
    ! coordinates to another.  The remapping can be chosen to be first-order 
    ! or second-order accurate.
    !
    ! Note: The cost of this subroutine scales as nlyr; a previous version scaled as nlyr^2.
    !
    ! Author: William Lipscomb, LANL

    use glide_types, only: HO_VERTICAL_REMAP_SECOND_ORDER

    implicit none
 
    ! in-out arguments
 
    integer, intent(in) ::  &
         nx, ny,     &! number of cells in EW and NS directions
         nlyr,       &! number of vertical layers
         ntracer      ! number of tracer fields

    type(parallel_type), intent(in) :: &
         parallel    ! info for parallel communication

    integer, intent(in) ::  &
         itest, jtest, rtest    ! coordinates of diagnostic point

    real(dp), dimension (nlyr+1), intent(in) ::  &
         sigma        ! sigma vertical coordinate (at layer interfaces)

    real(dp), dimension (nx, ny, nlyr), intent(inout) ::  &
         hlyr         ! layer thickness

    real(dp), dimension (nx, ny, ntracer, nlyr), intent(inout) ::   &
         trcr         ! tracer field to be remapped
                      ! tracer(k) = value at midpoint of layer k

    real(dp), dimension (nx, ny, ntracer), intent(in), optional ::   &
         trcr_usrf,  &! tracer field at upper surface
         trcr_lsrf    ! tracer field at lower surface

    integer, intent(in) ::  &
         vert_remap_accuracy    ! order of accuracy for vertical remapping
                                ! HO_VERTICAL_REMAP_FIRST_ORDER or HO_VERTICAL_REMAP_SECOMD_ORDER

    ! local variables
 
    integer :: i, j, k, k1, k2
    integer :: iglobal, jglobal
    integer :: nt

    real(dp), dimension(nlyr+1) ::  &
         z1,        &! layer interfaces in old coordinate system
                     ! z1(1) = 0. = value at upper surface
                     ! z1(k) = value at top of layer k
                     ! z1(nlyr+1) = value at lower surface (= 1 in sigma coordinates)
         z2          ! layer interfaces in new coordinate system

    real(dp), dimension(0:nlyr+1) ::  &
         zmid        ! layer midpoints in old coordinate system
                     ! zmid(0) = value at upper surface
                     ! zmid(nlyr+1) = value at lower surface

    real(dp) ::        &
         thck,      &! total thickness
         rthck       ! reciprocal of total thickness
 
    real(dp), dimension(ntracer,nlyr) ::       &
         gradt,     &! gradient of a tracer within a layer
         htsum       ! sum of thickness*tracer in a layer         

    real(dp), dimension(ntracer) ::  &
         tm1, tp1,       &! tracer values in layer k-1 and k+1
         tmax, tmin,     &! max/min tracer value in a layer and its neighbors
         tzmax, tzmin,   &! max/min value of reconstructed tracer in a layer (relative to midpt value)
         wk1, wk2         ! work arrays

    real(dp) :: zlo, zhi, dz, zav, hovlp

    character(len=256) :: message

    !WHL - debug
    real(dp) :: local_max, global_max
    real(dp) :: local_min, global_min

    !Note: Currently, temperature or enthalpy has tracer index 1; ice age, if computed, has index 2.
    !      Could make this more robust by having an nt_ice_age variable assigned at startup.
    integer, parameter :: nt_ice_age = 2

    ! Global max and min ice age
    if (verbose_ice_age) then
       local_min = minval(trcr(:,:,nt_ice_age,:))
       global_min = parallel_reduce_min(local_min)
       local_max = maxval(trcr(:,:,nt_ice_age,:))
       global_max = parallel_reduce_max(local_max)
       if (this_rank == rtest) then
          print*, 'Vertical remap, accuracy =', vert_remap_accuracy
          print*, 'max, min(ice_age):', global_max, global_min
       endif
    endif

    ! loop over locally owned cells
    do j = 1+nhalo, ny-nhalo
       do i = 1+nhalo, nx-nhalo

          !-----------------------------------------------------------------
          ! Compute total thickness
          !-----------------------------------------------------------------
          
          thck = 0.d0
          do k = 1, nlyr
             thck = thck + hlyr(i,j,k)
          enddo
          
          !-----------------------------------------------------------------
          ! If thck > 0, do vertical remapping of tracers
          ! WHL - Use tiny instead of 0.0d0 to avoid a rare divzero error.
          !-----------------------------------------------------------------

          if (thck > tiny(0.0d0)) then

             !-----------------------------------------------------------------
             ! Determine vertical coordinate z1, given input layer thicknesses.
             ! These are the coordinates from which we start.
             !-----------------------------------------------------------------

             rthck = 1.d0/thck

             z1(1) = 0.d0
             do k = 2, nlyr
                z1(k) = z1(k-1) +  hlyr(i,j,k-1)*rthck 
             enddo
             z1(nlyr+1) = 1.d0
                       
             !-----------------------------------------------------------------
             ! Compute layer midpoints for z1. Tracers are located at midpoints.
             !-----------------------------------------------------------------

             zmid(0) = 0.0d0
             do k = 1, nlyr
                zmid(k) = 0.5d0 * (z1(k) + z1(k+1))
             enddo
             zmid(nlyr+1) = z1(nlyr+1)
                
             !-----------------------------------------------------------------
             ! Compute vertical coordinate z2, given sigma.
             ! These are the coordinates to which we remap in the vertical.
             !-----------------------------------------------------------------
             
             z2(1) = 0.d0
             do k = 2, nlyr
                z2(k) = sigma(k)
             enddo
             z2(nlyr+1) = 1.d0

             if (verbose_ice_age .and. i == itest .and. j == jtest .and. this_rank == rtest) then
                write(6,*) ' '
                write(6,*) 'Column age, r, i, j =', rtest, i, j
                write(6,*) 'k, z1, dz, init ice_age (yr):'
                do k = 1,nlyr
                   write(6,'(i4,3f14.8)') k, z1(k), hlyr(i,j,k), trcr(i,j,nt_ice_age,k)*tim0/scyr
                enddo
             endif

             !-----------------------------------------------------------------
             ! Compute new layer thicknesses (z2 coordinates)
             !-----------------------------------------------------------------
             
             do k = 1, nlyr
                hlyr(i,j,k) = (z2(k+1) - z2(k)) * thck
             enddo
             
             !-----------------------------------------------------------------
             ! For second-order remapping: Compute tracer gradients in each layer.
             !                             Limit the gradients to preserve monotonicity.
             ! For first-order remapping:  Set the gradients to zero.
             !-----------------------------------------------------------------

             if (vert_remap_accuracy == HO_VERTICAL_REMAP_SECOND_ORDER) then
             
                do k = 1, nlyr
                   
                   ! Load values in layers above and below
                   
                   if (k > 1) then
                      tm1(:) = trcr(i,j,:,k-1)
                   else
                      if (present(trcr_usrf)) then
                         tm1(:) = trcr_usrf(i,j,:)
                      else
                         tm1(:) = trcr(i,j,:,1)
                      endif
                   endif
                   
                   if (k < nlyr) then
                      tp1(:) = trcr(i,j,:,k+1)
                   else
                      if (present(trcr_lsrf)) then
                         tp1(:) = trcr_lsrf(i,j,:)
                      else
                         tp1(:) = trcr(i,j,:,nlyr)
                      endif
                   endif
                   
                   ! Compute unlimited gradient
                   
                   dz = zmid(k+1) - zmid(k-1)
                   if (dz > 0.0d0) then
                      gradt(:,k) = (tp1(:) - tm1(:)) / dz
                   else
                      gradt(:,k) = 0.0d0
                   endif

                   ! Find the max and min deviations of tracer values in layer k-1, k, k+1
                   ! from the value in layer k
                   tmax(:) = max(tm1(:),trcr(i,j,:,k),tp1(:)) - trcr(i,j,:,k)
                   tmin(:) = min(tm1(:),trcr(i,j,:,k),tp1(:)) - trcr(i,j,:,k)
                
                   ! Find the max and min tracer deviations in this layer, given the unlimited gradient
                   wk1(:) = gradt(:,k) * (z1(k+1) - zmid(k))
                   wk2(:) = gradt(:,k) * (z1(k) - zmid(k))
                   tzmax(:) = max(wk1(:), wk2(:))
                   tzmin(:) = min(wk1(:), wk2(:))

                   ! Limit the gradient
                   
                   where (abs(tzmin) > 0.0d0)
                      wk1 = max(0.0d0, tmin/tzmin)
                   elsewhere
                      wk1 = 1.0d0
                   endwhere
                   
                   where (abs(tzmax) > 0.0d0)
                      wk2 = max(0.0d0, tmax/tzmax)
                   elsewhere
                      wk2 = 1.0d0
                   endwhere
                   
                   gradt(:,k) = gradt(:,k) * min(1.0d0, wk1(:), wk2(:))
                
                enddo   ! k
                
             else   ! first-order vertical remapping
                
                gradt(:,:) = 0.0d0
                
             endif

             !-----------------------------------------------------------------
             ! Compute sum of h*T for each new layer (k2) by integrating
             ! over the regions of overlap with old layers (k1).
             !
             ! The basic formula is as follows:
             !
             ! int_zlo^zhi [T(z) dz] = int_zlo^zhi [Tmid + gradT*(z - zmid)] dz
             !                       = dz * [Tmid + gradT*(zav - zmid)]
             ! where dz = zhi - zlo, zav = (zhi+zlo)/2, Tmid = midpoint tracer value
             !
             ! For first-order remapping, gradT = 0.
             !-----------------------------------------------------------------
          
             htsum(:,:) = 0.d0
             k1 = 1
             k2 = 1
             do while (k1 <= nlyr .and. k2 <= nlyr)
                zhi = min (z1(k1+1), z2(k2+1))
                zlo = max (z1(k1),   z2(k2))
                zav = 0.5d0 * (zlo + zhi)
                hovlp = max (zhi-zlo, 0.d0) * thck
                htsum(:,k2) = htsum(:,k2)   &
                            + hovlp * (trcr(i,j,:,k1) + gradt(:,k1) * (zav - zmid(k1)))
                if (z1(k1+1) > z2(k2+1)) then
                   k2 = k2 + 1
                else
                   k1 = k1 + 1
                endif
             enddo
          
             ! compute new tracer values
             ! WHL - Use tiny instead of 0.0d0 to avoid a rare divzero error.
             !       Such an error occurred during a 2-km Antarctic spin-up in Dec. 2019.
             ! Note: Since thck > tiny, we should have hlyr > tiny for all k.
             !       To be safe, however, allow for the possibility that thck > tiny but hlyr is not.

             do k = 1, nlyr
                if (hlyr(i,j,k) > tiny(0.0d0)) then
                   trcr(i,j,:,k) = htsum(:,k) / hlyr(i,j,k)
                else
                   trcr(i,j,:,k) = 0.0d0
                endif
             enddo

          else   ! thck = 0.0

             trcr(i,j,:,:) = 0.d0

          endif  ! thck > 0
             
       enddo       ! i
    enddo          ! j

    !WHL - debug - check for NaNs
    !TODO - Remove this check when satisfied the 2nd order remapping is working
    do k = 1, nlyr
       do nt = 1, ntracer
          do j = 1, ny
             do i = 1, nx
                if (trcr(i,j,nt,k) /= trcr(i,j,nt,k)) then
                   call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                   write(message,*) 'ERROR: Vertical remap, iglobal, jglobal, k, hlyr, trcr:', &
                        iglobal, jglobal, k, hlyr(i,j,k), trcr(i,j,nt,k)
                   call write_log(trim(message), GM_FATAL) 
                endif
             enddo
          enddo
       enddo
    enddo

    if (verbose_ice_age .and. this_rank == rtest .and. ntracer >= 2) then
       ! print column diagnostics
       i = itest
       j = jtest
       write(6,*) ' '
       write(6,*) 'k, z2, hlyr, remapped ice_age (yr):'
       do k = 1, nlyr
          write(6,'(i4,3f14.8)') k, z2(k), hlyr(i,j,k), trcr(i,j,nt_ice_age,k)*tim0/scyr
       enddo
    endif

    end subroutine glissade_vertical_remap

!=======================================================================


    subroutine upwind_field (nx,       ny,         &
                             ilo, ihi, jlo, jhi,   &
                             dx,       dy,         &
                             dt,       phi,        &
                             uee,      vnn)
      !
      ! first-order upwind transport algorithm
      !
      !
      ! Authors: Elizabeth Hunke and William Lipscomb, LANL
      !
      ! input/output arguments

      integer, intent (in) ::     &
         nx, ny             ,&! block dimensions
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real(dp), intent(in) ::         &
         dx, dy             ,&! x and y gridcell dimensions
         dt                   ! time step

      real(dp), dimension(nx,ny), &
         intent(inout) ::                       &
         phi                  ! scalar field

      real(dp), dimension(nx,ny),         &
         intent(in)::     &
         uee, vnn             ! cell edge velocities

      ! local variables

      integer ::     &
         i, j                   ! standard indices

      real(dp) ::        &
         upwind, y1, y2, a, h   ! function

      real(dp), dimension(nx,ny) ::  &
         worka, workb

    !-------------------------------------------------------------------
    ! Define upwind function
    !-------------------------------------------------------------------

      upwind(y1,y2,a,h) = 0.5d0*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

    !-------------------------------------------------------------------
    ! upwind transport
    !-------------------------------------------------------------------

      do j = jlo-1, jhi
      do i = ilo-1, ihi
         worka(i,j) = upwind(phi(i,j),phi(i+1,j),uee(i,j),dy)
         workb(i,j) = upwind(phi(i,j),phi(i,j+1),vnn(i,j),dx)
      enddo
      enddo

      do j = jlo, jhi
      do i = ilo, ihi
         phi(i,j) = phi(i,j) - ( worka(i,j)-worka(i-1,j)      &
                               + workb(i,j)-workb(i,j-1) )    &
                               / (dx*dy)
      enddo
      enddo

    end subroutine upwind_field

!=======================================================================

  end module glissade_transport

!=======================================================================
