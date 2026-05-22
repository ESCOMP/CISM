!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_diagnostics.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains routines for computing various diagnostic quantities.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  module glissade_diagnostics

    ! This module contains several subroutines for computing diagnostic quantities.
    ! Some are called at the end of each timestep, while others are called for specific MIPs,
    !  such as calvingMIP.
    ! TODO: Move the contents of glide_diagnostics to this module.

    use glimmer_global, only: dp
    use glimmer_paramets, only: iulog, eps11
    use glimmer_physcon, only: scyr
    !!    use glimmer_log
    use glimmer_utils, only: point_diag
    use glide_types
    use cism_parallel, only: this_rank, main_task, nhalo

    implicit none

    private
    public :: glissade_mass_balance_diagnostics, glissade_grounding_line_flux, &
         glissade_stress_tensor_eigenvalues, glissade_strain_rate_tensor_eigenvalues, &
         glissade_calvingmip_diag

    logical, parameter :: verbose_calvingmip = .true.

  contains

!****************************************************************************

  subroutine glissade_mass_balance_diagnostics(model)

    ! Compute mass balance diagnostics associated with five processes:
    ! (1) surface mass balance, (2) basal mass balance, (3) calving,
    ! (4) lateral melt, and (5) other kinds of ice removal.
    ! Ice cap removal, if applied, falls under (5).

    use glimmer_physcon, only: rhoi, rhow

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! Compute diagnostics

    ! surface mass balance in units of mm/yr w.e.
    ! (model%climate%acab has units of m/s of ice)
    ! Note: This is not necessary (and can break exact restart) if the SMB was already input in units of mm/yr
    if (model%options%smb_input /= SMB_INPUT_MMYR_WE) then
       model%climate%smb(:,:) = (model%climate%acab(:,:) * scyr) * (1000.d0 * rhoi/rhow)
    endif

    ! Corrections for basal melt at the calving front; convert basal melt to calving in CF cells.
    ! Computed melt rates can be large in CF cells when applying a calving mask and adjusting deltaT_ocn
    !  based on a thickness target.  In this case, it is better to think of the melt as part of the calving.
    ! Note: Both calving_thck and bmlt_applied have dimensionless model units;
    !       calving_thck = calving thickness per timestep, while bmlt_applied = melt per unit time

    if (model%options%whichcalving == CALVING_GRID_MASK .or. model%options%apply_calving_mask) then
       where (model%calving%calving_front_mask == 1)
          model%calving%calving_thck = model%calving%calving_thck + model%basal_melt%bmlt_applied * model%numerics%dt
          model%basal_melt%bmlt_applied = 0.0d0
       endwhere
    endif

    ! surface, basal, calving, lateral melt, and ice removal mass fluxes (kg/m^2/s)
    ! positive for mass gain, negative for mass loss
    model%mass_flux%sfc_mbal_flux(:,:) = rhoi * model%climate%acab_applied(:,:)
    model%mass_flux%basal_mbal_flux(:,:) = rhoi * (-model%basal_melt%bmlt_applied(:,:))
    model%mass_flux%calving_flux(:,:) = rhoi * (-model%calving%calving_thck(:,:)) / model%numerics%dt
    model%mass_flux%latmelt_flux(:,:) = rhoi * (-model%lateral_melt%melt_thck(:,:)) / model%numerics%dt
    model%mass_flux%removal_flux(:,:) = rhoi * (-model%geometry%removal_thck(:,:)) / model%numerics%dt

    ! rates of calving, lateral melt and ice removal (m/yr ice; positive for ice loss)
    model%calving%calving_rate(:,:) = model%calving%calving_thck(:,:) / (model%numerics%dt/scyr)
    model%lateral_melt%melt_rate(:,:) = model%lateral_melt%melt_thck(:,:) / (model%numerics%dt/scyr)
    model%geometry%removal_rate(:,:) = model%geometry%removal_thck(:,:) / (model%numerics%dt/scyr)

  end subroutine glissade_mass_balance_diagnostics

!---------------------------------------------------------------------------

  subroutine glissade_grounding_line_flux(&
       nx,                       ny,            &
       dx,                       dy,            &
       sigma,                                   &
       thck,                                    &
       uvel,                     vvel,          &
       ice_mask,                 floating_mask, &
       ocean_mask,                              &
       gl_flux_east,             gl_flux_north, &
       gl_flux)

    ! Compute northward and eastward land ice fluxes at grounding lines,
    !  and a cell-based grounding-line flux field.
    ! Note: Since the GL thicknesses are approximated, the GL fluxes will not exactly
    !        match the fluxes computed by the transport scheme.
    !       Also, the GL fluxes do not include thinning/calving of grounded marine cliffs.

    implicit none

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::                     &
         nx, ny                                   !> horizontal grid dimensions

    real(dp), intent(in) ::                    &
         dx, dy                                   !> horizontal grid spacing

    real(dp), dimension(:), intent(in) ::      &
         sigma                                    !> vertical sigma coordinate

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck                                     !> ice thickness

    real(dp), dimension(:,:,:), intent(in) ::  &
         uvel, vvel                               !> ice velocity in x and y directions

    integer, dimension(nx,ny), intent(in) ::  &
         ice_mask,                              & !> = 1 where ice is present, else = 0
         floating_mask,                         & !> = 1 where ice is present and floating, else = 0
         ocean_mask                               !> = 1 for ice-free ocean, else = 0

    ! Note: gl_flux_east and gl_flux_north are directional
    !       (positive for eastward/northward, negative for westward/southward)
    !       gl_flux is a cell-based quantity based on flux magnitudes on each edge
    !       (so gl_flux >= 0)

    real(dp), dimension(:,:), intent(out) ::   &
         gl_flux_east,                          & !> grounding line flux on east edges
         gl_flux_north,                         & !> grounding line flux on north edges
         gl_flux                                  !> grounding line flux per grid cell

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer  :: i,j,k                                     !> local cell indices
    integer  :: upn                                       !> vertical grid dimension
    real(dp), dimension(:), allocatable :: uavg, vavg     !> local horizontal velocity averages
    real(dp) :: thck_gl                                   !> GL thickness derived from topg_gl

    upn = size(sigma)

    allocate(uavg(upn), vavg(upn))

    ! Initialize
    gl_flux_east(:,:)  = 0.d0
    gl_flux_north(:,:) = 0.d0
    gl_flux(:,:)       = 0.d0

    ! Compute grounding line fluxes on east and north edges.
    ! Look for edges with a grounded cell on one side and a floating cell on the other.

    do j = nhalo+1, ny-nhalo
        do i = nhalo+1, nx-nhalo

            ! check east edge
           if ( (   (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) .and.   &  ! (i,j) grounded
                (ocean_mask(i+1,j) == 1 .or.  floating_mask(i+1,j) == 1) )     &  ! (i+1,j) floating or ocean
                                        .or.                                   &
                ( (ice_mask(i+1,j) == 1 .and. floating_mask(i+1,j) == 0) .and. &  ! (i+1,j) grounded
                  (ocean_mask(i,j) == 1  .or. floating_mask(i,j) == 1) ) ) then   ! (i,j) floating or ocean

                uavg(:) = (uvel(:,i,j) + uvel(:,i,j-1)) / 2.d0
                if (ice_mask(i,j) == 1 .and. ice_mask(i+1,j) == 1) then
                   ! set GL thickness to the average thickness of the two cells
                   thck_gl = (thck(i,j) + thck(i+1,j)) / 2.d0
                else
                   ! set GL thickness to the thickness of the ice-filled cell
                   thck_gl = max(thck(i,j), thck(i+1,j))
                endif

                do k = 1, upn-1
                    gl_flux_east(i,j) = gl_flux_east(i,j) &
                                        + thck_gl * (sigma(k+1) - sigma(k)) * (uavg(k) + uavg(k+1))/2.d0
                enddo
            endif

            ! check north edge
           if ( (   (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) .and.   &  ! (i,j) grounded
                (ocean_mask(i,j+1) == 1 .or.  floating_mask(i,j+1) == 1) )     &  ! (i,j+1) floating or ocean
                                        .or.                                   &
                ( (ice_mask(i,j+1) == 1 .and. floating_mask(i,j+1) == 0) .and. &  ! (i,j+1) grounded
                  (ocean_mask(i,j) == 1  .or. floating_mask(i,j) == 1) ) ) then   ! (i,j) floating or ocean

                vavg(:) = (vvel(:,i-1,j) + vvel(:,i,j)) / 2.d0
                if (ice_mask(i,j) == 1 .and. ice_mask(i,j+1) == 1) then
                   ! set GL thickness to the average thickness of the two cells
                   thck_gl = (thck(i,j) + thck(i,j+1)) / 2.d0
                else
                   ! set GL thickness to the thickness of the ice-filled cell
                   thck_gl = max(thck(i,j), thck(i,j+1))
                endif

                do k = 1, upn-1
                    gl_flux_north(i,j) = gl_flux_north(i,j) &
                                        + thck_gl * (sigma(k+1) - sigma(k)) * (vavg(k) + vavg(k+1))/2.d0
                enddo
             endif

        enddo   ! i
    enddo   ! j

    ! Compute mass flux through grounding line in each cell.
    ! Only a grounded cell can lose mass. We need to check the direction of the fluxes.

    do j = nhalo+1,ny-nhalo
        do i = nhalo+1,nx-nhalo

            ! Check the sign for east-west flow and assign the flux accordingly
            if (gl_flux_east(i,j) < 0.d0) then
                ! The ice is flowing westward and the flux belongs to the right adjacent cell
                gl_flux(i+1,j) = gl_flux(i+1,j) - gl_flux_east(i,j)
            else
                ! The ice is flowing eastward and the flux belongs to this cell
                gl_flux(i,j) = gl_flux(i,j) + gl_flux_east(i,j)
            endif

            ! Check the sign for north-south flow and assign the flux accordingly
            if (gl_flux_north(i,j) < 0.d0) then
                ! The ice is flowing southward and the flux belongs to the top adjacent cell
                gl_flux(i,j+1) = gl_flux(i,j+1) - gl_flux_north(i,j)
            else
                ! The ice is flowing northward and the flux belongs to this cell
                gl_flux(i,j) = gl_flux(i,j) + gl_flux_north(i,j)
            endif

        enddo   ! i
    enddo   ! j

    ! Convert from m^2/s to kg/m/s
    gl_flux_east  = gl_flux_east  * rhoi
    gl_flux_north = gl_flux_north * rhoi
    gl_flux       = gl_flux       * rhoi

    deallocate(uavg, vavg)

  end subroutine glissade_grounding_line_flux

!---------------------------------------------------------------------------

  subroutine glissade_stress_tensor_eigenvalues(&
       nx,    ny,   nz,   &
       sigma,             &
       tau,               &
       tau_eigen1,        &
       tau_eigen2)

    ! Diagnose the eigenvalues of the 2D horizontal stress tensor.
    ! These are used for eigencalving and damage-based calving.

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
             tau_xx = tau_xx + tau%xx(k,i,j) * dsigma
             tau_yy = tau_yy + tau%yy(k,i,j) * dsigma
             tau_xy = tau_xy + tau%xy(k,i,j) * dsigma
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

    ! Diagnose the eigenvalues of the 2D horizontal strain rate tensor.
    ! These can be used for eigencalving and damage-based calving, or for diagnostics.
    ! There are two ways to call the subroutine:
    ! (1) Pass in the strain rate tensor and compute the eigenvalues directly.
    ! (2) Pass in the stress tensor as an optional argument, compute the strain rate tensor
    !     from the stress tensor and effective viscosity, and then compute the eigenvalues.

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
          strain_rate%scalar = tau%scalar / (2.d0 * efvs)
          strain_rate%xz = tau%xz / (2.d0 * efvs)
          strain_rate%yz = tau%yz / (2.d0 * efvs)
          strain_rate%xx = tau%xx / (2.d0 * efvs)
          strain_rate%yy = tau%yy / (2.d0 * efvs)
          strain_rate%xy = tau%xy / (2.d0 * efvs)
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
! The next three subroutines are diagnostic subroutines for CalvingMIP.
! They estimate the calving front location along 8 prescribed axes
!  for the circular and Thule domains.
! They are not necessary if we have offline tools for locating the CF,
!  but are included for flexibility.
! For details, see the CalvingMIP Wiki: https://github.com/JRowanJordan/CalvingMIP/wiki    
! See also the paper by J. Jordan et al. (2026, TC).
!---------------------------------------------------------------------------

  subroutine glissade_calvingmip_diag(model)

    ! Compute diagnostics for the CalvingMIP experiments
    ! These include the calvingMIP location, ice speed, and ice thickness
    !  along 8 axes for the circular and Thule domains.

    use glissade_masks, only: glissade_get_masks, glissade_calving_front_mask
    use glissade_grid_operators, only: glissade_unstagger
    use cism_parallel, only: parallel_halo, parallel_global_sum

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! Local variables

    integer :: n
    integer :: nx, ny
    integer :: itest, jtest, rtest
    real(dp) :: dx, dy

    type(parallel_type) :: parallel   ! info for parallel communication

    integer, dimension(model%general%ewn, model%general%nsn) :: &
         ice_mask,                & ! = 1 if ice is present
         floating_mask,           & ! = 1 if ice is present and floating
         land_mask,               & ! = 1 if topg - eus >= 0
         ocean_mask                 ! = 1 if ice is absent and topg - eus < 0

    integer, dimension(model%general%ewn, model%general%nsn) :: &
         partial_cf_mask,         & ! = 1 for partially filled CF cells (thck < thck_effective), else = 0
         full_mask                  ! = 1 for ice-filled cells that are not partial_cf cells, else = 0

    real(dp), dimension(model%general%ewn-1,model%general%nsn-1) :: &
         velnorm_mean               ! mean ice speed at vertices (m/s)

    integer, dimension(model%general%ewn-1,model%general%nsn-1) :: &
         vmask                      ! = 1 for vertices of active cells

    real(dp), dimension(model%general%ewn,model%general%nsn) :: &
         uvel,                    & ! uvel_2d averaged to cell centers (m/s)
         vvel                       ! vvel_2d averaged to cell centers (m/s)

    real(dp) :: &
         total_ice_area         ! total effective ice area (with weighting by effective_areafrac)

    real(dp), dimension(4) :: quadrant_sum   ! sum over each of the 4 quadrants for calvingMIP

    nx = model%general%ewn
    ny = model%general%nsn

    dx = model%numerics%dew
    dy = model%numerics%dns

    parallel = model%parallel

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    ! Compute the ice speed at cell centers, averaged from neighboring vertices.
    ! Include in the average only vertices with nonzero speeds (i.e., ice present)

    velnorm_mean = sqrt(model%velocity%uvel_2d**2 + model%velocity%vvel_2d**2)

    where (velnorm_mean > 0.0d0)
       vmask = 1
    elsewhere
       vmask = 0
    endwhere

    ! Interpolate the velocity from cell vertices to centers.
    ! 'stagger_margin_in = 1' means that masked-out values are not part of the average.

    call glissade_unstagger(&
         nx,                     ny,      &
         model%velocity%uvel_2d, uvel,    &
         vmask,                  stagger_margin_in = 1)

    call glissade_unstagger(&
         nx,                     ny,      &
         model%velocity%vvel_2d, vvel,    &
         vmask,                  stagger_margin_in = 1)

    call parallel_halo(uvel, parallel)
    call parallel_halo(vvel, parallel)

    ! Compute some required masks

    call glissade_get_masks(&
         nx,            ny,             &
         parallel,                      &
         model%geometry%thck,           &
         model%geometry%topg,           &
         model%climate%eus,             &
         eps11,                         &
         ice_mask,                      &
         floating_mask = floating_mask, &
         ocean_mask = ocean_mask,       &
         land_mask = land_mask)

    call glissade_calving_front_mask(&
         nx,            ny,                &
         model%options%which_ho_calving_front,   &
         parallel,                         &
         itest,  jtest, rtest,             &
         model%geometry%thck,              &
         model%geometry%topg,              &
         model%climate%eus,                &
         ice_mask,      floating_mask,     &
         ocean_mask,    land_mask,         &
         model%calving%calving_front_mask, &
         model%calving%dthck_dx_cf,        &
         dx,            dy,                &
         model%calving%thck_effective,     &
         model%calving%thck_effective_min, &
         partial_cf_mask,                  &
         full_mask,                        &
         model%calving%effective_areafrac)

    ! Compute the diagnostics for the chosen domain (circular or Thule)
    if (model%options%which_ho_calvingmip_domain == HO_CALVINGMIP_DOMAIN_CIRCULAR .and. &
        model%options%whichcalving == CF_ADVANCE_RETREAT_RATE) then

       call locate_calving_front_circular(&
            nx,             ny,                &
            dx,             dy,                &  ! m
            model%general%x0,                  &  ! m
            model%general%y0,                  &  ! m
            model%general%x1,                  &  ! m
            model%general%y1,                  &  ! m
            parallel,                          &
            itest,   jtest,   rtest,           &
            model%calving%effective_areafrac,  &
            model%calving%thck_effective,      &
            uvel,           vvel,              &  ! m/s
            model%calving%cf_locx,             &  ! m
            model%calving%cf_locy,             &  ! m
            model%calving%cf_radius,           &  ! m
            model%calving%cf_thck,             &  ! m
            model%calving%cf_uvel,             &  ! m/s
            model%calving%cf_vvel)

    elseif (model%options%which_ho_calvingmip_domain == HO_CALVINGMIP_DOMAIN_THULE .and. &
            model%options%whichcalving == CF_ADVANCE_RETREAT_RATE) then

       call locate_calving_front_thule(&
            nx,             ny,                &
            dx,             dy,                &  ! m
            model%general%x0,                  &  ! m
            model%general%y0,                  &  ! m
            model%general%x1,                  &  ! m
            model%general%y1,                  &  ! m
            parallel,                          &
            itest,   jtest,   rtest,           &
            model%calving%effective_areafrac,  &
            model%calving%thck_effective,      &
            uvel,           vvel,              &  ! m/s
            model%calving%cf_locx,             &  ! m
            model%calving%cf_locy,             &  ! m
            model%calving%cf_radius,           &  ! m
            model%calving%cf_thck,             &  ! m
            model%calving%cf_uvel,             &  ! m/s
            model%calving%cf_vvel)

    endif

    if (verbose_calvingmip) then

       ! Compute the total ice area and the area of each quadrant
       total_ice_area = parallel_global_sum(dx*dy*model%calving%effective_areafrac, parallel)

       call sum_over_quadrants(&
            nx,                  ny,                &
            parallel,                               &
            model%calving%effective_areafrac, &  ! m^2
            quadrant_sum)

       if (this_rank == rtest) then
          write(iulog,*) 'Total ice area (km^2):', total_ice_area/1.0d6
          write(iulog,*) 'Quadrant area:'
          do n = 1, 4
             write(iulog,*) n, quadrant_sum(n)
          enddo
       endif

    endif

  end subroutine glissade_calvingmip_diag

!---------------------------------------------------------------------------

  subroutine locate_calving_front_circular(&
       nx,               ny,           &
       dx,               dy,           &
       x0,               y0,           &
       x1,               y1,           &
       parallel,                       &
       itest,   jtest,   rtest,        &
       areafrac,                       &
       thck_effective,                 &
       uvel,             vvel,         &
       cf_locx,          cf_locy,      &
       cf_radius,        cf_thck,      &
       cf_uvel,          cf_vvel)

    use cism_parallel, only: parallel_reduce_maxloc, parallel_reduce_minloc, broadcast
    use glissade_grid_operators, only: glissade_stagger

    ! Find the calving front location along eight profiles on the circular domain.
    ! These profiles are the four cardinal directions (N, S, E, W) along with the diagonals
    !  that form 45-degree angles with the cardinal directions.
    !
    ! Method for finding the CF location along the x or y axis:
    ! (1) Identify the last cell (i,j) with areafrac >= 0.5, followed by the first cell with areafrac < 0.5.
    ! (2) Interpolate linearly between the two cells to find the point where areafrac = 0.5.
    ! (3) Do a similar interpolation for CF thickness (using thck_effective) and velocity components.
    !
    ! The method for diagonal axes is similar, except that we add cell corners to the interpolation.
    ! The CF is located either (1) between a cell center with areafrac >= 0.5 and a corner with areafrac < 0.5,
    ! or (2) between a corner with areafrac >= 0.5 and the following center with areafrac < 0.5.

    integer, intent(in) :: &
         nx, ny,                 & ! grid dimensions
         itest, jtest, rtest       ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dx, dy                    ! grid cell size (m)

    real(dp), dimension(nx-1), intent(in) :: x0  ! x coordinate of NE cell corners
    real(dp), dimension(ny-1), intent(in) :: y0  ! y coordinate of NE cell corners
    real(dp), dimension(nx), intent(in) :: x1    ! x coordinate of cell centers
    real(dp), dimension(ny), intent(in) :: y1    ! y coordinate of cell centers

    type(parallel_type), intent(in) :: &
         parallel                    ! info for parallel communication

    real(dp), dimension(nx,ny), intent(in) :: &
         areafrac,                 & ! effective fractional area, in range [0,1]
         thck_effective,           & ! ice thickness (m)
         uvel, vvel                  ! ice velocity components (m/s)

    real(dp), dimension(8), intent(out) :: &
         cf_locx, cf_locy,         & ! x and y components of CF location (m)
         cf_radius,                & ! radial distance of CF from origin (m)
         cf_thck,                  & ! ice thickness at CF (m)
         cf_uvel, cf_vvel            ! u and v velocity components at CF (m/s)

    ! local variables

    integer :: i, j, iglobal, jglobal
    integer :: axis
    integer :: procnum
    real(dp) :: cf_loc_xmax, cf_loc_ymax, cf_loc_xmin, cf_loc_ymin
    real(dp) :: wt_factor
    real(dp) :: this_areafrac, next_areafrac, corner_frac
    real(dp) :: this_thck, next_thck, corner_thck
    real(dp) :: this_uvel, next_uvel, corner_uvel
    real(dp) :: this_vvel, next_vvel, corner_vvel
    real(dp) :: speed

    ! Note: This subroutine assumes that the grid origin (0,0) is located at a cell vertex,
    !       so the x- and y-axes lie along cell edges.

    logical :: &
         x_axis_thru_edges,    & ! true if the x-axis passes through cell edges
         y_axis_thru_edges       ! true if the y-axis passes through cell edges

    ! Find the x and y coordinates of the calving front along the different axes
    !  specified in CalvingMIP.

    if (this_rank == rtest) write(iulog,*) 'Locate_calving_front for calvingMIP, rtest =', rtest

    ! Determine whether the x and/or y axes pass through cell edges on this processor
    x_axis_thru_edges = .false.
    do j = nhalo+1, ny-nhalo
       if (y0(j) == 0.0d0) then
          x_axis_thru_edges = .true.
          exit
       endif
    enddo

    y_axis_thru_edges = .false.
    do i = nhalo+1, nx-nhalo
       if (x0(i) == 0.0d0) then
          y_axis_thru_edges = .true.
          exit
       endif
    enddo

    if (verbose_calvingmip) then
       if (x_axis_thru_edges) then
!          write(iulog,*) this_rank, 'x_axis_thru_edges', x_axis_thru_edges
       endif
       if (y_axis_thru_edges) then
!          write(iulog,*) this_rank, 'y_axis_thru_edges', y_axis_thru_edges
       endif
    endif

    ! Initialize calvingMIP diagnostics
    cf_locx = 0.0d0
    cf_locy = 0.0d0
    cf_radius = 0.0d0
    cf_thck = 0.0d0
    cf_uvel = 0.0d0
    cf_vvel = 0.0d0

    ! Find the CF location along each of 8 axes
    ! The CF has areafrac = 0.5. To find its location, interpolate linearly between
    !  two points, one with areafrac >= 0.5 and one with areafrac < 0.5.
    ! All loops are over locally owned cells

    ! Compute diagnostics for profiles along the x- and y-axes

    axis = 1  ! index for the positive y-axis (profile A)
    if (y_axis_thru_edges) then
       do i = nhalo+1, nx-nhalo
          if (x0(i) == 0.0d0) then  ! E edge of cell lies on the y-axis
             cf_locx(axis) = 0.0d0
             do j = nhalo+1, ny-nhalo
                this_areafrac = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
                next_areafrac = 0.5d0 * (areafrac(i,j+1) + areafrac(i+1,j+1))
                if (this_areafrac >= 0.5d0 .and. next_areafrac < 0.5d0) then
                   ! CF lies between j and j+1
                   wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                   cf_locy(axis) = y1(j)*wt_factor + y1(j+1)*(1.0d0 - wt_factor)
                   cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                   if (next_areafrac > eps11) then  ! take a weighted average between j and j+1
                      this_thck = 0.5d0 * (thck_effective(i,j) + thck_effective(i+1,j))
                      next_thck = 0.5d0 * (thck_effective(i,j+1) + thck_effective(i+1,j+1))
                      cf_thck(axis) = this_thck*wt_factor + next_thck*(1.0d0 - wt_factor)
                      this_uvel = 0.5d0 * (uvel(i,j) + uvel(i+1,j))
                      next_uvel = 0.5d0 * (uvel(i,j+1) + uvel(i+1,j+1))
                      cf_uvel(axis) = this_uvel*wt_factor + next_uvel*(1.0d0 - wt_factor)
                      this_vvel = 0.5d0 * (vvel(i,j) + vvel(i+1,j))
                      next_vvel = 0.5d0 * (vvel(i,j+1) + vvel(i+1,j+1))
                      cf_vvel(axis) = this_vvel*wt_factor + next_vvel*(1.0d0 - wt_factor)
                   else  ! use the values from this j
                      cf_thck(axis) = 0.5d0 * (thck_effective(i,j) + thck_effective(i+1,j))
                      cf_uvel(axis) = 0.5d0 * (uvel(i,j) + uvel(i+1,j))
                      cf_vvel(axis) = 0.5d0 * (vvel(i,j) + vvel(i+1,j))
                   endif
                endif
             enddo
          endif
       enddo
    endif   ! y_axis_thru_edges

    ! If this proc has a positive value of y, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_locy(axis), xout=cf_loc_ymax, xprocout=procnum)

    ! Broadcast the calvingMIP axis 1 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)


    axis = 3  ! index for the positive x-axis (profile C)
    if (x_axis_thru_edges) then
       do j = nhalo+1, ny-nhalo
          if (y0(j) == 0.0d0) then
             cf_locy(axis) = 0.0d0
             do i = nhalo+1, nx-nhalo
                this_areafrac = 0.5d0 * (areafrac(i,j) + areafrac(i,j+1))
                next_areafrac = 0.5d0 * (areafrac(i+1,j) + areafrac(i+1,j+1))
                if (this_areafrac >= 0.5d0 .and. next_areafrac < 0.5d0) then
                   ! CF lies between i and i+1
                   wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                   cf_locx(axis) = x1(i)*wt_factor + x1(i+1)*(1.0d0 - wt_factor)
                   cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                   if (next_areafrac > eps11) then  ! take a weighted average between i and i+1
                      this_thck = 0.5d0 * (thck_effective(i,j) + thck_effective(i,j+1))
                      next_thck = 0.5d0 * (thck_effective(i+1,j) + thck_effective(i+1,j+1))
                      cf_thck(axis) = this_thck*wt_factor + next_thck*(1.0d0 - wt_factor)
                      this_uvel = 0.5d0 * (uvel(i,j) + uvel(i,j+1))
                      next_uvel = 0.5d0 * (uvel(i+1,j) + uvel(i+1,j+1))
                      cf_uvel(axis) = this_uvel*wt_factor + next_uvel*(1.0d0 - wt_factor)
                      this_vvel = 0.5d0 * (vvel(i,j) + vvel(i,j+1))
                      next_vvel = 0.5d0 * (vvel(i+1,j) + vvel(i+1,j+1))
                      cf_vvel(axis) = this_vvel*wt_factor + next_vvel*(1.0d0 - wt_factor)
                   else  ! use the values from this i
                      cf_thck(axis) = 0.5d0 * (thck_effective(i,j) + thck_effective(i,j+1))
                      cf_uvel(axis) = 0.5d0 * (uvel(i,j) + uvel(i,j+1))
                      cf_vvel(axis) = 0.5d0 * (vvel(i,j) + vvel(i,j+1))
                   endif
                endif
             enddo
          endif
       enddo
    endif   ! x_axis_thru_edges


    ! If this proc has a positive value of x, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_locx(axis), xout=cf_loc_xmax, xprocout=procnum)

    ! Broadcast the calvingMIP axis 3 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)

    axis = 5  ! index for the negative y-axis (profile E)
    if (y_axis_thru_edges) then
       do i = nhalo+1, nx-nhalo
          if (x0(i) == 0.0d0) then  ! E edge of cell lies on the y-axis
             cf_locx(axis) = 0.0d0
             do j = ny-nhalo, nhalo+1, -1
                this_areafrac = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
                next_areafrac = 0.5d0 * (areafrac(i,j-1) + areafrac(i+1,j-1))
                if (this_areafrac >= 0.5d0 .and. next_areafrac < 0.5d0) then
                   ! CF lies between j and j-1
                   wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                   cf_locy(axis) = y1(j)*wt_factor + y1(j-1)*(1.0d0 - wt_factor)
                   cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                   if (next_areafrac > eps11) then  ! take a weighted average between j and j-1
                      this_thck = 0.5d0 * (thck_effective(i,j) + thck_effective(i+1,j))
                      next_thck = 0.5d0 * (thck_effective(i,j-1) + thck_effective(i+1,j-1))
                      cf_thck(axis) = this_thck*wt_factor + next_thck*(1.0d0 - wt_factor)
                      this_uvel = 0.5d0 * (uvel(i,j) + uvel(i+1,j))
                      next_uvel = 0.5d0 * (uvel(i,j-1) + uvel(i+1,j-1))
                      cf_uvel(axis) = this_uvel*wt_factor + next_uvel*(1.0d0 - wt_factor)
                      this_vvel = 0.5d0 * (vvel(i,j) + vvel(i+1,j))
                      next_vvel = 0.5d0 * (vvel(i,j-1) + vvel(i+1,j-1))
                      cf_vvel(axis) = this_vvel*wt_factor + next_vvel*(1.0d0 - wt_factor)
                   else  ! use the values from this j
                      cf_thck(axis) = 0.5d0 * (thck_effective(i,j) + thck_effective(i+1,j))
                      cf_uvel(axis) = 0.5d0 * (uvel(i,j) + uvel(i+1,j))
                      cf_vvel(axis) = 0.5d0 * (vvel(i,j) + vvel(i+1,j))
                   endif
                endif
             enddo
          endif
       enddo
    endif   ! y_axis_thru_edges

    ! If this proc has a negative value of y, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_locy(axis), xout=cf_loc_ymin, xprocout=procnum)

    ! Broadcast the calvingMIP axis 5 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)


    axis = 7  ! index for the negative x-axis (profile G)
    if (x_axis_thru_edges) then
       do j = nhalo+1, ny-nhalo
          if (y0(j) == 0.0d0) then
             cf_locy(axis) = 0.0d0
             do i = nx-nhalo, nhalo+1, -1
                this_areafrac = 0.5d0 * (areafrac(i,j) + areafrac(i,j+1))
                next_areafrac = 0.5d0 * (areafrac(i-1,j) + areafrac(i-1,j+1))
                if (this_areafrac >= 0.5d0 .and. next_areafrac < 0.5d0) then
                   ! CF lies between i and i-1
                   wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                   cf_locx(axis) = x1(i)*wt_factor + x1(i-1)*(1.0d0 - wt_factor)
                   cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                   if (next_areafrac > eps11) then  ! take a weighted average between i and i+1
                      this_thck = 0.5d0 * (thck_effective(i,j) + thck_effective(i,j+1))
                      next_thck = 0.5d0 * (thck_effective(i-1,j) + thck_effective(i-1,j+1))
                      cf_thck(axis) = this_thck*wt_factor + next_thck*(1.0d0 - wt_factor)
                      this_uvel = 0.5d0 * (uvel(i,j) + uvel(i,j+1))
                      next_uvel = 0.5d0 * (uvel(i-1,j) + uvel(i-1,j+1))
                      cf_uvel(axis) = this_uvel*wt_factor + next_uvel*(1.0d0 - wt_factor)
                      this_vvel = 0.5d0 * (vvel(i,j) + vvel(i,j+1))
                      next_vvel = 0.5d0 * (vvel(i-1,j) + vvel(i-1,j+1))
                      cf_vvel(axis) = this_vvel*wt_factor + next_vvel*(1.0d0 - wt_factor)
                   else  ! use the values from this i
                      cf_thck(axis) = 0.5d0 * (thck_effective(i,j) + thck_effective(i,j+1))
                      cf_uvel(axis) = 0.5d0 * (uvel(i,j) + uvel(i,j+1))
                      cf_vvel(axis) = 0.5d0 * (vvel(i,j) + vvel(i,j+1))
                   endif
                endif
             enddo
          endif
       enddo
    endif   ! x_axis_thru_edges

    ! If this proc has a negative value of x, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_locx(axis), xout=cf_loc_xmin, xprocout=procnum)

    ! Broadcast the calvingMIP axis 7 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)


    ! Compute diagnostics for the four profiles along the diagonals (y = x and y = -x)

    axis = 2  ! index for the line y = x in the positive x and y direction (profile B)
    do i = nhalo+1, nx-nhalo
       do j = nhalo+1, ny-nhalo
          if (x1(i) == y1(j)) then ! on the line y = x
             corner_frac = 0.5d0*(areafrac(i,j+1) + areafrac(i+1,j))
             if (areafrac(i,j) >= 0.5d0 .and. corner_frac < 0.5d0) then
                ! CF lies in the NE quadrant of cell (i,j)
                this_areafrac = areafrac(i,j)
                next_areafrac = corner_frac
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locx(axis) = x1(i)*wt_factor + x0(i)*(1.0d0 - wt_factor)
                cf_locy(axis) = y1(j)*wt_factor + y0(j)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                if (corner_frac > eps11) then  ! take a weighted average of the center and corner values
                   corner_thck = 0.5d0 * (thck_effective(i,j+1) + thck_effective(i+1,j))
                   cf_thck(axis) = thck_effective(i,j)*wt_factor + corner_thck*(1.0d0 - wt_factor)
                   corner_uvel = 0.5d0 * (uvel(i,j+1) + uvel(i+1,j))
                   cf_uvel(axis) = uvel(i,j)*wt_factor + corner_uvel*(1.0d0 - wt_factor)
                   corner_vvel = 0.5d0 * (vvel(i,j+1) + vvel(i+1,j))
                   cf_vvel(axis) = vvel(i,j)*wt_factor + corner_vvel*(1.0d0 - wt_factor)
                else  ! use the center values
                   cf_thck(axis) = thck_effective(i,j)
                   cf_uvel(axis) = uvel(i,j)
                   cf_vvel(axis) = vvel(i,j)
                endif
             elseif (corner_frac >= 0.5d0 .and. areafrac(i+1,j+1) < 0.5d0) then
                ! CF lies in the SW quadrant of cell (i+1,j+1)
                this_areafrac = corner_frac
                next_areafrac = areafrac(i+1,j+1)
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locx(axis) = x0(i)*wt_factor + x1(i+1)*(1.0d0 - wt_factor)
                cf_locy(axis) = y0(j)*wt_factor + y1(j+1)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                if (areafrac(i+1,j+1) > eps11) then  ! take a weighted average of the corner and center values
                   corner_thck = 0.5d0 * (thck_effective(i,j+1) + thck_effective(i+1,j))
                   cf_thck(axis) = corner_thck*wt_factor + thck_effective(i+1,j+1)*(1.0d0 - wt_factor)
                   corner_uvel = 0.5d0 * (uvel(i,j+1) + uvel(i+1,j))
                   cf_uvel(axis) = corner_uvel*wt_factor + uvel(i+1,j+1)*(1.0d0 - wt_factor)
                   corner_vvel = 0.5d0 * (vvel(i,j+1) + vvel(i+1,j))
                   cf_vvel(axis) = corner_vvel*wt_factor + vvel(i+1,j+1)*(1.0d0 - wt_factor)
                else   ! use the corner values
                   cf_thck(axis) = 0.5d0 * (thck_effective(i,j+1) + thck_effective(i+1,j))
                   cf_uvel(axis) = 0.5d0 * (uvel(i,j+1) + uvel(i+1,j))
                   cf_vvel(axis) = 0.5d0 * (vvel(i,j+1) + vvel(i+1,j))
                endif
             endif
          endif   ! on the line y = x
       enddo   ! i
    enddo   ! j

    ! If this proc has a positive value of x, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_locx(axis), xout=cf_loc_xmax, xprocout=procnum)

    ! Broadcast the calvingMIP axis 2 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)


    axis = 4  ! index for the line y = -x in the positive x and negative y direction (profile D)
    do i = nhalo+1, nx-nhalo
       do j = nhalo+1, ny-nhalo
          if (x1(i) == (-1.0d0)*y1(j)) then ! on the line y = -x
             corner_frac = 0.5d0*(areafrac(i,j-1) + areafrac(i+1,j))
             if (areafrac(i,j) >= 0.5d0 .and. corner_frac < 0.5d0) then
                ! CF lies in the SE quadrant of cell (i,j)
                this_areafrac = areafrac(i,j)
                next_areafrac = corner_frac
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locx(axis) = x1(i)*wt_factor + x0(i)*(1.0d0 - wt_factor)
                cf_locy(axis) = y1(j)*wt_factor + y0(j-1)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                if (corner_frac > eps11) then  ! take a weighted average of the center and corner values
                   corner_thck = 0.5d0 * (thck_effective(i,j-1) + thck_effective(i+1,j))
                   cf_thck(axis) = thck_effective(i,j)*wt_factor + corner_thck*(1.0d0 - wt_factor)
                   corner_uvel = 0.5d0 * (uvel(i,j-1) + uvel(i+1,j))
                   cf_uvel(axis) = uvel(i,j)*wt_factor + corner_uvel*(1.0d0 - wt_factor)
                   corner_vvel = 0.5d0 * (vvel(i,j-1) + vvel(i+1,j))
                   cf_vvel(axis) = vvel(i,j)*wt_factor + corner_vvel*(1.0d0 - wt_factor)
                else  ! use the center values
                   cf_thck(axis) = thck_effective(i,j)
                   cf_uvel(axis) = uvel(i,j)
                   cf_vvel(axis) = vvel(i,j)
                endif
             elseif (corner_frac >= 0.5d0 .and. areafrac(i+1,j-1) < 0.5d0) then
                ! CF lies in the NW quadrant of cell (i+1,j-1)
                this_areafrac = corner_frac
                next_areafrac = areafrac(i+1,j-1)
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locx(axis) = x0(i)*wt_factor + x1(i+1)*(1.0d0 - wt_factor)
                cf_locy(axis) = y0(j-1)*wt_factor + y1(j-1)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                if (areafrac(i+1,j-1) > eps11) then  ! take a weighted average of the corner and center values
                   corner_thck = 0.5d0 * (thck_effective(i,j-1) + thck_effective(i+1,j))
                   cf_thck(axis) = corner_thck*wt_factor + thck_effective(i+1,j-1)*(1.0d0 - wt_factor)
                   corner_uvel = 0.5d0 * (uvel(i,j-1) + uvel(i+1,j))
                   cf_uvel(axis) = corner_uvel*wt_factor + uvel(i+1,j-1)*(1.0d0 - wt_factor)
                   corner_vvel = 0.5d0 * (vvel(i,j-1) + vvel(i+1,j))
                   cf_vvel(axis) = corner_vvel*wt_factor + vvel(i+1,j-1)*(1.0d0 - wt_factor)
                else   ! use the corner values
                   cf_thck(axis) = 0.5d0 * (thck_effective(i,j-1) + thck_effective(i+1,j))
                   cf_uvel(axis) = 0.5d0 * (uvel(i,j-1) + uvel(i+1,j))
                   cf_vvel(axis) = 0.5d0 * (vvel(i,j-1) + vvel(i+1,j))
                endif
             endif
          endif
       enddo   ! i
    enddo   ! j

    ! If this proc has a positive value of x, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_locx(axis), xout=cf_loc_xmax, xprocout=procnum)

    ! Broadcast the calvingMIP axis 4 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)


    axis = 6  ! index for the line y = x in the negative x and y direction (profile F)
    do i = nhalo+1, nx-nhalo
       do j = nhalo+1, ny-nhalo
          if (x1(i) == y1(j)) then ! on the line y = x
             corner_frac = 0.5d0*(areafrac(i,j-1) + areafrac(i-1,j))
             if (areafrac(i,j) >= 0.5d0 .and. corner_frac < 0.5d0) then
                ! CF lies in the SW quadrant of cell (i,j)
                this_areafrac = areafrac(i,j)
                next_areafrac = corner_frac
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locx(axis) = x1(i)*wt_factor + x0(i-1)*(1.0d0 - wt_factor)
                cf_locy(axis) = y1(j)*wt_factor + y0(j-1)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                if (corner_frac > eps11) then  ! take a weighted average of the center and corner values
                   corner_thck = 0.5d0 * (thck_effective(i,j-1) + thck_effective(i-1,j))
                   cf_thck(axis) = thck_effective(i,j)*wt_factor + corner_thck*(1.0d0 - wt_factor)
                   corner_uvel = 0.5d0 * (uvel(i,j-1) + uvel(i-1,j))
                   cf_uvel(axis) = uvel(i,j)*wt_factor + corner_uvel*(1.0d0 - wt_factor)
                   corner_vvel = 0.5d0 * (vvel(i,j-1) + vvel(i-1,j))
                   cf_vvel(axis) = vvel(i,j)*wt_factor + corner_vvel*(1.0d0 - wt_factor)
                else  ! use the center values
                   cf_thck(axis) = thck_effective(i,j)
                   cf_uvel(axis) = uvel(i,j)
                   cf_vvel(axis) = vvel(i,j)
                endif
             elseif (corner_frac >= 0.5d0 .and. areafrac(i-1,j-1) < 0.5d0) then
                ! CF lies in the NE quadrant of cell (i-1,j-1)
                this_areafrac = corner_frac
                next_areafrac = areafrac(i-1,j-1)
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locx(axis) = x0(i-1)*wt_factor + x1(i-1)*(1.0d0 - wt_factor)
                cf_locy(axis) = y0(j-1)*wt_factor + y1(j-1)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                if (areafrac(i-1,j-1) > eps11) then  ! take a weighted average of the corner and center values
                   corner_thck = 0.5d0 * (thck_effective(i,j-1) + thck_effective(i-1,j))
                   cf_thck(axis) = corner_thck*wt_factor + thck_effective(i-1,j-1)*(1.0d0 - wt_factor)
                   corner_uvel = 0.5d0 * (uvel(i,j-1) + uvel(i-1,j))
                   cf_uvel(axis) = corner_uvel*wt_factor + uvel(i-1,j-1)*(1.0d0 - wt_factor)
                   corner_vvel = 0.5d0 * (vvel(i,j-1) + vvel(i-1,j))
                   cf_vvel(axis) = corner_vvel*wt_factor + vvel(i-1,j-1)*(1.0d0 - wt_factor)
                else   ! use the corner values
                   cf_thck(axis) = 0.5d0 * (thck_effective(i,j-1) + thck_effective(i-1,j))
                   cf_uvel(axis) = 0.5d0 * (uvel(i,j-1) + uvel(i-1,j))
                   cf_vvel(axis) = 0.5d0 * (vvel(i,j-1) + vvel(i-1,j))
                endif
             endif
          endif
       enddo   ! i
    enddo   ! j

    ! If this proc has a negative value of x, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_locx(axis), xout=cf_loc_xmin, xprocout=procnum)

    ! Broadcast the calvingMIP axis 6 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)


    axis = 8  ! index for the line y = -x in the negative x and positive y direction (profile H)
    do i = nhalo+1, nx-nhalo
       do j = nhalo+1, ny-nhalo
          if (x1(i) == (-1.0d0)*y1(j)) then ! on the line y = -x
             corner_frac = 0.5d0*(areafrac(i,j+1) + areafrac(i-1,j))
             if (areafrac(i,j) >= 0.5d0 .and. corner_frac < 0.5d0) then
                ! CF lies in the NW quadrant of cell (i,j)
                this_areafrac = areafrac(i,j)
                next_areafrac = corner_frac
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locx(axis) = x1(i)*wt_factor + x0(i-1)*(1.0d0 - wt_factor)
                cf_locy(axis) = y1(j)*wt_factor + y0(j)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                if (corner_frac > eps11) then  ! take a weighted average of the center and corner values
                   corner_thck = 0.5d0 * (thck_effective(i,j+1) + thck_effective(i-1,j))
                   cf_thck(axis) = thck_effective(i,j)*wt_factor + corner_thck*(1.0d0 - wt_factor)
                   corner_uvel = 0.5d0 * (uvel(i,j+1) + uvel(i-1,j))
                   cf_uvel(axis) = uvel(i,j)*wt_factor + corner_uvel*(1.0d0 - wt_factor)
                   corner_vvel = 0.5d0 * (vvel(i,j+1) + vvel(i-1,j))
                   cf_vvel(axis) = vvel(i,j)*wt_factor + corner_vvel*(1.0d0 - wt_factor)
                else  ! use the center values
                   cf_thck(axis) = thck_effective(i,j)
                   cf_uvel(axis) = uvel(i,j)
                   cf_vvel(axis) = vvel(i,j)
                endif
             elseif (corner_frac >= 0.5d0 .and. areafrac(i-1,j+1) < 0.5d0) then
                ! CF lies in the SE quadrant of cell (i-1,j+1)
                this_areafrac = corner_frac
                next_areafrac = areafrac(i-1,j+1)
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locx(axis) = x0(i-1)*wt_factor + x1(i-1)*(1.0d0 - wt_factor)
                cf_locy(axis) = y0(j)*wt_factor + y1(j+1)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                if (areafrac(i-1,j+1) > eps11) then  ! take a weighted average of the corner and center values
                   corner_thck = 0.5d0 * (thck_effective(i,j+1) + thck_effective(i-1,j))
                   cf_thck(axis) = corner_thck*wt_factor + thck_effective(i-1,j+1)*(1.0d0 - wt_factor)
                   corner_uvel = 0.5d0 * (uvel(i,j+1) + uvel(i-1,j))
                   cf_uvel(axis) = corner_uvel*wt_factor + uvel(i-1,j+1)*(1.0d0 - wt_factor)
                   corner_vvel = 0.5d0 * (vvel(i,j+1) + vvel(i-1,j))
                   cf_vvel(axis) = corner_vvel*wt_factor + vvel(i-1,j+1)*(1.0d0 - wt_factor)
                else   ! use the corner values
                   cf_thck(axis) = 0.5d0 * (thck_effective(i,j+1) + thck_effective(i-1,j))
                   cf_uvel(axis) = 0.5d0 * (uvel(i,j+1) + uvel(i-1,j))
                   cf_vvel(axis) = 0.5d0 * (vvel(i,j+1) + vvel(i-1,j))
                endif
             endif
          endif
       enddo   ! i
    enddo   ! j

    ! If this proc has a negative value of x, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_locx(axis), xout=cf_loc_xmin, xprocout=procnum)

    ! Broadcast the calvingMIP axis 8 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)

    if (verbose_calvingmip .and. main_task) then
       write(iulog,*) ' '
       write(iulog,*) 'Circular domain: axis, x_cf, y_cf, radius (km), thck (m), uvel, vvel, speed (m/yr)'
       do axis = 1, 8
          speed = sqrt(cf_uvel(axis)**2 + cf_vvel(axis)**2)
          write(iulog,'(i4,7f15.8)') axis, cf_locx(axis)/1000.d0, cf_locy(axis)/1000.d0, &
               cf_radius(axis)/1000.d0, cf_thck(axis), cf_uvel(axis)*scyr, cf_vvel(axis)*scyr, speed*scyr
       enddo
    endif

  end subroutine locate_calving_front_circular

!---------------------------------------------------------------------------


  subroutine locate_calving_front_thule(&
       nx,               ny,           &
       dx,               dy,           &
       x0,               y0,           &
       x1,               y1,           &
       parallel,                       &
       itest,   jtest,   rtest,        &
       areafrac,                       &
       thck_effective,                 &
       uvel,             vvel,         &
       cf_locx,          cf_locy,      &
       cf_radius,        cf_thck,      &
       cf_uvel,          cf_vvel)

    use cism_parallel, only: parallel_reduce_maxloc, parallel_reduce_minloc, &
         broadcast, parallel_globalindex

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
         parallel                     ! info for parallel communication

    real(dp), dimension(nx,ny), intent(in) :: &
         areafrac,                  & ! effective fractional area, in range [0,1]
         thck_effective,            & ! effective ice thickness (m)
         uvel, vvel                   ! ice velocity components (m/s)

    ! Note: Axis 1 is Caprona A and axis 2 is Halbrane A; both are in the upper left (NW) quadrant
    real(dp), dimension(8), intent(out) :: &
         cf_locx, cf_locy,          & ! x and y components of CF location (m)
         cf_radius,                 & ! radial distance of CF from origin (m)
         cf_thck,                   & ! ice thickness at CF (m)
         cf_uvel, cf_vvel             ! u and v velocity components at CF (m/s)

    ! local variables

    integer :: i, j, jj, iglobal, jglobal
    integer :: axis
    integer :: procnum
    real(dp) :: cf_loc_xmax, cf_loc_ymax, cf_loc_xmin, cf_loc_ymin
    real(dp) :: wt_factor, wt1, wt2
    real(dp) :: this_areafrac, next_areafrac
    real(dp) :: this_thck, next_thck
    real(dp) :: this_uvel, next_uvel
    real(dp) :: this_vvel, next_vvel
    real(dp) :: x_intercept, y_intercept, slope  ! properties of the profile
    real(dp) :: x_lim, y_lim                     ! outer limits of the profile
    real(dp) :: speed

    real(dp), dimension(ny) ::  &
         areafrac_yint,          & ! areafrac at the point (x_int(j), y1(j))
         x_int                     ! x value of profile where it intersects with y1(j)

    ! Find the x and y coordinates of the calving front for the Thule domain
    ! along the different profiles specified in CalvingMIP.

    ! Initialize the output arrays
    cf_locx = 0.0d0
    cf_locy = 0.0d0
    cf_radius = 0.0d0
    cf_thck = 0.0d0
    cf_uvel = 0.0d0
    cf_vvel = 0.0d0

    ! Find the CF location along the different axes.
    ! The Caprona axes are labeled 1, 3, 5 and 7; Halbrane axes are 2, 4, 6 and 8.
    ! The axes labeled 1 and 2 are Caprona A and Halbrane A, as shown in the Jordon et al. paper
    ! All loops are over locally owned cells.
    ! Assume that the x and y axes coincide with cell edges (not cell centers).

    ! Compute diagnostics for the four Caprona profiles
    ! Note: The Caprona profiles cut across cells without passing through centers or corners.
    !       As a result, the logic below is more complicated than for the Halbrane profiles.

    axis = 1  ! index for Caprona A
    x_intercept = -390.d3
    x_lim = -590.d3
    y_lim = 450.d3
    slope = y_lim/(x_lim - x_intercept)  ! rise over run = 450/(-200) = -2.25
    y_intercept = -x_intercept * slope

    ! Adjust x_lim and y_lim to allow the CF to be a little out of bounds
    x_lim = x_lim * 1.2d0
    y_lim = y_lim * 1.2d0

    x_int = 0.0d0
    areafrac_yint = 0.0d0

    do j = nhalo, ny-nhalo+1
       if (y1(j) > 0.0d0 .and. y1(j) <= y_lim) then  ! y1 in range
          x_int(j) = (y1(j) - y_intercept)/slope  ! profile intersects y1 grid at (x_int(j),y1(j))
          do i = nx, 2, -1
             if (x_int(j) <= x1(i) .and. x_int(j) > x1(i-1)) then
                ! Interpolate to estimate areafrac at (x_int(j),y1(j))
                wt_factor = (x1(i) - x_int(j))/dx
                areafrac_yint(j) = (1.0d0 - wt_factor)*areafrac(i,j) + wt_factor*areafrac(i-1,j)
                exit
             endif
          enddo
       endif
    enddo

    do j = nhalo, ny-nhalo
       if (areafrac_yint(j) >= 0.5d0 .and. areafrac_yint(j+1) < 0.5d0) then
          ! Find the point along the axis where the interpolated areafrac = 0.5
          wt_factor = (areafrac_yint(j) - 0.5d0) / (areafrac_yint(j) - areafrac_yint(j+1))
          cf_locx(axis) = (1.0d0 - wt_factor)*x_int(j) + wt_factor*x_int(j+1)
          cf_locy(axis) = (1.0d0 - wt_factor)*y1(j) + wt_factor*y1(j+1)
          cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)

          ! Estimate thck and other variables at point (x_int(j), y1(j))
          do i = nx, 2, -1
             if (x_int(j) <= x1(i) .and. x_int(j) > x1(i-1)) then
                wt1 = (x1(i) - x_int(j))/dx
                if (areafrac(i,j) > eps11 .and. areafrac(i-1,j) > eps11) then
                   this_thck = (1.0d0 - wt1)*thck_effective(i,j) + wt1*thck_effective(i-1,j)
                   this_uvel = (1.0d0 - wt1)*uvel(i,j) + wt1*uvel(i-1,j)
                   this_vvel = (1.0d0 - wt1)*vvel(i,j) + wt1*vvel(i-1,j)
                elseif (areafrac(i,j) > eps11) then
                   this_thck = thck_effective(i,j)
                   this_uvel = uvel(i,j)
                   this_vvel = vvel(i,j)
                elseif (areafrac(i-1,j) > eps11) then
                   this_thck = thck_effective(i-1,j)
                   this_uvel = uvel(i-1,j)
                   this_vvel = vvel(i-1,j)
                else   ! this should not happen
                   this_thck = 0.0d0
                   this_uvel = 0.0d0
                   this_vvel = 0.0d0
                endif
                exit
             endif
          enddo   ! i

          ! Estimate thck and other variables at point (x_int(j+1), y1(j+1))
          do i = nx, 2, -1
             if (x_int(j+1) <= x1(i) .and. x_int(j+1) > x1(i-1)) then
                wt2 = (x1(i) - x_int(j+1))/dx
                if (areafrac(i,j+1) > eps11 .and. areafrac(i-1,j+1) > eps11) then
                   next_thck = (1.0d0 - wt2)*thck_effective(i,j+1) + wt2*thck_effective(i-1,j+1)
                   next_uvel = (1.0d0 - wt2)*uvel(i,j+1) + wt2*uvel(i-1,j+1)
                   next_vvel = (1.0d0 - wt2)*vvel(i,j+1) + wt2*vvel(i-1,j+1)
                elseif (areafrac(i,j+1) > eps11) then
                   next_thck = thck_effective(i,j+1)
                   next_uvel = uvel(i,j+1)
                   next_vvel = vvel(i,j+1)
                elseif (areafrac(i-1,j+1) > eps11) then
                   next_thck = thck_effective(i-1,j+1)
                   next_uvel = uvel(i-1,j+1)
                   next_vvel = vvel(i-1,j+1)
                else
                   next_thck = 0.0d0
                   next_uvel = 0.0d0
                   next_vvel = 0.0d0
                endif
                exit
             endif
          enddo   ! i

          !Note: Should have nonzero velocity wherever thck > 0
          if (this_thck > eps11 .and. next_thck > eps11) then
             cf_thck(axis) = (1.0d0 - wt_factor)*this_thck + wt_factor*next_thck
             cf_uvel(axis) = (1.0d0 - wt_factor)*this_uvel + wt_factor*next_uvel
             cf_vvel(axis) = (1.0d0 - wt_factor)*this_vvel + wt_factor*next_vvel
          elseif (this_thck > eps11) then
             cf_thck(axis) = this_thck
             cf_uvel(axis) = this_uvel
             cf_vvel(axis) = this_vvel
          elseif (next_thck > eps11) then
             cf_thck(axis) = next_thck
             cf_uvel(axis) = next_uvel
             cf_vvel(axis) = next_vvel
          else
             cf_thck(axis) = 0.0d0
             cf_uvel(axis) = 0.0d0
             cf_vvel(axis) = 0.0d0
          endif
          call parallel_globalindex(i, j, iglobal, jglobal, parallel)
!          write(iulog,*) 'Axis 1, possible CF: rank, i, j, ig, jg:', this_rank, i, j, iglobal, jglobal
!          write(iulog,*) '     y1(j), x_int(j), x_int(j+1), areafrac_int(j), areafrac_int(j+1):',  &
!               y1(j), x_int(j), x_int(j+1), areafrac_yint(j), areafrac_yint(j+1)
!          write(iulog,*) '     x_cf, y_cf:', cf_locx(axis), cf_locy(axis)
!          write(iulog,*) 'this_thck, next_thck, thck_cf:', this_thck, next_thck, cf_thck(axis)
       endif  ! areafrac_yint(j) >= 0.5, areafrac_yint(j+1) < 0.5
    enddo   ! j

    ! If this proc has a positive value of y, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_locy(axis), xout=cf_loc_ymax, xprocout=procnum)

    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)


    axis = 3  ! index for Caprona B
    x_intercept = 390.d3
    x_lim = 590.d3
    y_lim = 450.d3
    slope = y_lim/(x_lim - x_intercept)  ! rise over run = 450/200 = 2.25
    y_intercept = -x_intercept * slope

    ! Adjust x_lim and y_lim to allow the CF to be a little out of bounds
    x_lim = x_lim * 1.2d0
    y_lim = y_lim * 1.2d0

    x_int = 0.0d0
    areafrac_yint = 0.0d0

    do j = nhalo, ny-nhalo+1
       if (y1(j) > 0.0d0 .and. y1(j) <= y_lim) then  ! y1 in range
          x_int(j) = (y1(j) - y_intercept)/slope  ! profile intersects y1 grid at (x_int(j),y1(j))
          do i = 1, nx-1
             if (x_int(j) >= x1(i) .and. x_int(j) < x1(i+1)) then
                ! Interpolate to estimate areafrac at (x_int(j),y1(j))
                wt_factor = (x_int(j) - x1(i))/dx
                areafrac_yint(j) = (1.0d0 - wt_factor)*areafrac(i,j) + wt_factor*areafrac(i+1,j)
                exit
             endif
          enddo
       endif
    enddo

    do j = nhalo, ny-nhalo
       if (areafrac_yint(j) >= 0.5d0 .and. areafrac_yint(j+1) < 0.5d0) then
          ! Find the point along the axis where the interpolated areafrac = 0.5
          wt_factor = (areafrac_yint(j) - 0.5d0) / (areafrac_yint(j) - areafrac_yint(j+1))
          cf_locx(axis) = (1.0d0 - wt_factor)*x_int(j) + wt_factor*x_int(j+1)
          cf_locy(axis) = (1.0d0 - wt_factor)*y1(j) + wt_factor*y1(j+1)
          cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)

          ! Estimate thck and other variables at point (x_int(j), y1(j))
          do i = 1, nx-1
             if (x_int(j) >= x1(i) .and. x_int(j) < x1(i+1)) then
                wt1 = (x_int(j) - x1(i))/dx
                if (areafrac(i,j) > eps11 .and. areafrac(i+1,j) > eps11) then
                   this_thck = (1.0d0 - wt1)*thck_effective(i,j) + wt1*thck_effective(i+1,j)
                   this_uvel = (1.0d0 - wt1)*uvel(i,j) + wt1*uvel(i+1,j)
                   this_vvel = (1.0d0 - wt1)*vvel(i,j) + wt1*vvel(i+1,j)
                elseif (areafrac(i,j) > eps11) then
                   this_thck = thck_effective(i,j)
                   this_uvel = uvel(i,j)
                   this_vvel = vvel(i,j)
                elseif (areafrac(i+1,j) > eps11) then
                   this_thck = thck_effective(i+1,j)
                   this_uvel = uvel(i+1,j)
                   this_vvel = vvel(i+1,j)
                else   ! this should not happen
                   this_thck = 0.0d0
                   this_uvel = 0.0d0
                   this_vvel = 0.0d0
                endif
                exit
             endif
          enddo   ! i

          ! Estimate thck and other variables at point (x_int(j+1), y1(j+1))
          do i = 1, nx-1
             if (x_int(j+1) >= x1(i) .and. x_int(j+1) < x1(i+1)) then
                wt2 = (x_int(j+1) - x1(i))/dx
                if (areafrac(i,j+1) > eps11 .and. areafrac(i+1,j+1) > eps11) then
                   next_thck = (1.0d0 - wt2)*thck_effective(i,j+1) + wt2*thck_effective(i+1,j+1)
                   next_uvel = (1.0d0 - wt2)*uvel(i,j+1) + wt2*uvel(i+1,j+1)
                   next_vvel = (1.0d0 - wt2)*vvel(i,j+1) + wt2*vvel(i+1,j+1)
                elseif (areafrac(i,j+1) > eps11) then
                   next_thck = thck_effective(i,j+1)
                   next_uvel = uvel(i,j+1)
                   next_vvel = vvel(i,j+1)
                elseif (areafrac(i+1,j+1) > eps11) then
                   next_thck = thck_effective(i+1,j+1)
                   next_uvel = uvel(i+1,j+1)
                   next_vvel = vvel(i+1,j+1)
                else
                   next_thck = 0.0d0
                   next_uvel = 0.0d0
                   next_vvel = 0.0d0
                endif
                exit
             endif
          enddo   ! i

          !Note: Should have nonzero velocity wherever thck > 0
          if (this_thck > eps11 .and. next_thck > eps11) then
             cf_thck(axis) = (1.0d0 - wt_factor)*this_thck + wt_factor*next_thck
             cf_uvel(axis) = (1.0d0 - wt_factor)*this_uvel + wt_factor*next_uvel
             cf_vvel(axis) = (1.0d0 - wt_factor)*this_vvel + wt_factor*next_vvel
          elseif (this_thck > eps11) then
             cf_thck(axis) = this_thck
             cf_uvel(axis) = this_uvel
             cf_vvel(axis) = this_vvel
          elseif (next_thck > eps11) then
             cf_thck(axis) = next_thck
             cf_uvel(axis) = next_uvel
             cf_vvel(axis) = next_vvel
          else
             cf_thck(axis) = 0.0d0
             cf_uvel(axis) = 0.0d0
             cf_vvel(axis) = 0.0d0
          endif
       endif  ! areafrac_yint(j) >= 0.5, areafrac_yint(j+1) < 0.5
    enddo   ! j

    ! If this proc has a positive value of y, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_locy(axis), xout=cf_loc_ymax, xprocout=procnum)

    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)


    axis = 5  ! index for Caprona C
    x_intercept = -390.d3
    x_lim = -590.d3
    y_lim = -450.d3
    slope = y_lim/(x_lim - x_intercept)  ! rise over run = -450/(-200) = 9/4
    y_intercept = -x_intercept * slope

    ! Adjust x_lim and y_lim to allow the CF to be a little out of bounds
    x_lim = x_lim * 1.2d0
    y_lim = y_lim * 1.2d0

    x_int = 0.0d0
    areafrac_yint = 0.0d0

    do j = ny-nhalo+1, nhalo, -1
       if (y1(j) < 0.0d0 .and. y1(j) >= y_lim) then  ! y1 in range
          x_int(j) = (y1(j) - y_intercept)/slope  ! profile intersects y1 grid at (x_int(j),y1(j))
          do i = nx, 2, -1
             if (x_int(j) >= x1(i-1) .and.  x_int(j) < x1(i)) then
                ! Interpolate to estimate areafrac at (x_int(j),y1(j))
                wt_factor = (x1(i) - x_int(j))/dx
                areafrac_yint(j) = (1.0d0 - wt_factor)*areafrac(i,j) + wt_factor*areafrac(i-1,j)
                exit
             endif
          enddo
       endif
    enddo

    do j = ny-nhalo+1, nhalo+1, -1
       if (areafrac_yint(j) >= 0.5d0 .and. areafrac_yint(j-1) < 0.5d0) then
          ! Find the point along the axis where the interpolated areafrac = 0.5
          wt_factor = (areafrac_yint(j) - 0.5d0) / (areafrac_yint(j) - areafrac_yint(j-1))
          cf_locx(axis) = (1.0d0 - wt_factor)*x_int(j) + wt_factor*x_int(j-1)
          cf_locy(axis) = (1.0d0 - wt_factor)*y1(j) + wt_factor*y1(j-1)
          cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)

          ! Estimate thck and other variables at point (x_int(j), y1(j))
          do i = nx, 2, -1
             if (x_int(j) <= x1(i) .and. x_int(j) > x1(i-1)) then
                wt1 = (x1(i) - x_int(j))/dx
                if (areafrac(i,j) > eps11 .and. areafrac(i-1,j) > eps11) then
                   this_thck = (1.0d0 - wt1)*thck_effective(i,j) + wt1*thck_effective(i-1,j)
                   this_uvel = (1.0d0 - wt1)*uvel(i,j) + wt1*uvel(i-1,j)
                   this_vvel = (1.0d0 - wt1)*vvel(i,j) + wt1*vvel(i-1,j)
                elseif (areafrac(i,j) > eps11) then
                   this_thck = thck_effective(i,j)
                   this_uvel = uvel(i,j)
                   this_vvel = vvel(i,j)
                elseif (areafrac(i-1,j) > eps11) then
                   this_thck = thck_effective(i-1,j)
                   this_uvel = uvel(i-1,j)
                   this_vvel = vvel(i-1,j)
                else   ! this should not happen
                   this_thck = 0.0d0
                   this_uvel = 0.0d0
                   this_vvel = 0.0d0
                endif
                exit
             endif
          enddo   ! i

          ! Estimate thck and other variables at point (x_int(j+1), y1(j-1))
          do i = nx, 2, -1
             if (x_int(j-1) <= x1(i) .and. x_int(j-1) > x1(i-1)) then
                wt2 = (x1(i) - x_int(j-1))/dx
                if (areafrac(i,j-1) > eps11 .and. areafrac(i-1,j-1) > eps11) then
                   next_thck = (1.0d0 - wt2)*thck_effective(i,j-1) + wt2*thck_effective(i-1,j-1)
                   next_uvel = (1.0d0 - wt2)*uvel(i,j-1) + wt2*uvel(i-1,j-1)
                   next_vvel = (1.0d0 - wt2)*vvel(i,j-1) + wt2*vvel(i-1,j-1)
                elseif (areafrac(i,j-1) > eps11) then
                   next_thck = thck_effective(i,j-1)
                   next_uvel = uvel(i,j-1)
                   next_vvel = vvel(i,j-1)
                elseif (areafrac(i-1,j-1) > eps11) then
                   next_thck = thck_effective(i-1,j-1)
                   next_uvel = uvel(i-1,j-1)
                   next_vvel = vvel(i-1,j-1)
                else
                   next_thck = 0.0d0
                   next_uvel = 0.0d0
                   next_vvel = 0.0d0
                endif
                exit
             endif
          enddo   ! i

          !Note: Should have nonzero velocity wherever thck > 0
          if (this_thck > eps11 .and. next_thck > eps11) then
             cf_thck(axis) = (1.0d0 - wt_factor)*this_thck + wt_factor*next_thck
             cf_uvel(axis) = (1.0d0 - wt_factor)*this_uvel + wt_factor*next_uvel
             cf_vvel(axis) = (1.0d0 - wt_factor)*this_vvel + wt_factor*next_vvel
          elseif (this_thck > eps11) then
             cf_thck(axis) = this_thck
             cf_uvel(axis) = this_uvel
             cf_vvel(axis) = this_vvel
          elseif (next_thck > eps11) then
             cf_thck(axis) = next_thck
             cf_uvel(axis) = next_uvel
             cf_vvel(axis) = next_vvel
          else
             cf_thck(axis) = 0.0d0
             cf_uvel(axis) = 0.0d0
             cf_vvel(axis) = 0.0d0
          endif

       endif  ! areafrac_yint(j) >= 0.5, areafrac_yint(j-1) < 0.5
    enddo   ! j

    ! If this proc has a negative value of y, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_locy(axis), xout=cf_loc_ymin, xprocout=procnum)

    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)
       

    axis = 7  ! index for Caprona D
    x_intercept = 390.d3
    x_lim = 590.d3
    y_lim = -450.d3
    slope = y_lim/(x_lim - x_intercept)  ! rise over run = -450/200 = -9/4
    y_intercept = -x_intercept * slope

    ! Adjust x_lim and y_lim to allow the CF to be a little out of bounds
    x_lim = x_lim * 1.2d0
    y_lim = y_lim * 1.2d0

    x_int = 0.0d0
    areafrac_yint = 0.0d0

    do j = ny-nhalo+1, nhalo, -1
       if (y1(j) < 0.0d0 .and. y1(j) >= y_lim) then  ! y1 in range
          x_int(j) = (y1(j) - y_intercept)/slope  ! profile intersects y1 grid at (x_int(j),y1(j))
          do i = 1, nx-1
             if (x_int(j) >= x1(i) .and. x_int(j) < x1(i+1)) then
                ! Interpolate to estimate areafrac at (x_int(j),y1(j))
                wt_factor = (x_int(j) - x1(i))/dx
                areafrac_yint(j) = (1.0d0 - wt_factor)*areafrac(i,j) + wt_factor*areafrac(i+1,j)
                exit
             endif
          enddo
       endif
    enddo

    do j = ny-nhalo+1, nhalo+1, -1
       if (areafrac_yint(j) >= 0.5d0 .and. areafrac_yint(j-1) < 0.5d0) then
          ! Find the point along the axis where the interpolated areafrac = 0.5
          wt_factor = (areafrac_yint(j) - 0.5d0) / (areafrac_yint(j) - areafrac_yint(j-1))
          cf_locx(axis) = (1.0d0 - wt_factor)*x_int(j) + wt_factor*x_int(j-1)
          cf_locy(axis) = (1.0d0 - wt_factor)*y1(j) + wt_factor*y1(j-1)
          cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)

          ! Estimate thck and other variables at point (x_int(j), y1(j))
          do i = 1, nx-1
             if (x_int(j) >= x1(i) .and. x_int(j) < x1(i+1)) then
                wt1 = (x_int(j) - x1(i))/dx
                if (areafrac(i,j) > eps11 .and. areafrac(i+1,j) > eps11) then
                   this_thck = (1.0d0 - wt1)*thck_effective(i,j) + wt1*thck_effective(i+1,j)
                   this_uvel = (1.0d0 - wt1)*uvel(i,j) + wt1*uvel(i+1,j)
                   this_vvel = (1.0d0 - wt1)*vvel(i,j) + wt1*vvel(i+1,j)
                elseif (areafrac(i,j) > eps11) then
                   this_thck = thck_effective(i,j)
                   this_uvel = uvel(i,j)
                   this_vvel = vvel(i,j)
                elseif (areafrac(i+1,j) > eps11) then
                   this_thck = thck_effective(i+1,j)
                   this_uvel = uvel(i+1,j)
                   this_vvel = vvel(i+1,j)
                else   ! this should not happen
                   this_thck = 0.0d0
                   this_uvel = 0.0d0
                   this_vvel = 0.0d0
                endif
                exit
             endif
          enddo   ! i

          ! Estimate thck and other variables at point (x_int(j+1), y1(j-1))
          do i = 1, nx-1
             if (x_int(j-1) >= x1(i) .and. x_int(j-1) < x1(i+1)) then
                wt2 = (x_int(j-1) - x1(i))/dx
                if (areafrac(i,j-1) > eps11 .and. areafrac(i+1,j-1) > eps11) then
                   next_thck = (1.0d0 - wt2)*thck_effective(i,j-1) + wt2*thck_effective(i+1,j-1)
                   next_uvel = (1.0d0 - wt2)*uvel(i,j-1) + wt2*uvel(i+1,j-1)
                   next_vvel = (1.0d0 - wt2)*vvel(i,j-1) + wt2*vvel(i+1,j-1)
                elseif (areafrac(i,j-1) > eps11) then
                   next_thck = thck_effective(i,j-1)
                   next_uvel = uvel(i,j-1)
                   next_vvel = vvel(i,j-1)
                elseif (areafrac(i+1,j-1) > eps11) then
                   next_thck = thck_effective(i+1,j-1)
                   next_uvel = uvel(i+1,j-1)
                   next_vvel = vvel(i+1,j-1)
                else
                   next_thck = 0.0d0
                   next_uvel = 0.0d0
                   next_vvel = 0.0d0
                endif
                exit
             endif
          enddo   ! i

          !Note: Should have nonzero velocity wherever thck > 0
          if (this_thck > eps11 .and. next_thck > eps11) then
             cf_thck(axis) = (1.0d0 - wt_factor)*this_thck + wt_factor*next_thck
             cf_uvel(axis) = (1.0d0 - wt_factor)*this_uvel + wt_factor*next_uvel
             cf_vvel(axis) = (1.0d0 - wt_factor)*this_vvel + wt_factor*next_vvel
          elseif (this_thck > eps11) then
             cf_thck(axis) = this_thck
             cf_uvel(axis) = this_uvel
             cf_vvel(axis) = this_vvel
          elseif (next_thck > eps11) then
             cf_thck(axis) = next_thck
             cf_uvel(axis) = next_uvel
             cf_vvel(axis) = next_vvel
          else
             cf_thck(axis) = 0.0d0
             cf_uvel(axis) = 0.0d0
             cf_vvel(axis) = 0.0d0
          endif

       endif  ! areafrac_yint(j) >= 0.5, areafrac_yint(j-1) < 0.5
    enddo   ! j

    ! If this proc has a negative value of y, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_locy(axis), xout=cf_loc_ymin, xprocout=procnum)

    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)

    

    ! Compute diagnostics for the four Halbrane profiles
    ! Find a point along each profile where the interpolated areafrac = 0.5

    axis = 2  ! index for Halbrane A
    x_intercept = -150.d3

    do i = nhalo, nx-nhalo
       if (abs(x0(i) - x_intercept) < eps11) then  ! E edge of cell lies on the vertical Halbrane profile
          cf_locx(axis) = x_intercept
          do j = nhalo, ny-nhalo
             this_areafrac = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
             next_areafrac = 0.5d0 * (areafrac(i,j+1) + areafrac(i+1,j+1))
             if (this_areafrac >= 0.5d0 .and. next_areafrac < 0.5d0) then
                ! CF lies between j and j+1
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locy(axis) = y1(j)*wt_factor + y1(j+1)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                ! The following logic allows for the possibility that one of the two neighbor cells is ice-free
                if (areafrac(i,j) > eps11 .and. areafrac(i+1,j) > eps11) then
                   this_thck = 0.5d0 * (thck_effective(i,j) + thck_effective(i+1,j))
                   this_uvel = 0.5d0 * (uvel(i,j) + uvel(i+1,j))
                   this_vvel = 0.5d0 * (vvel(i,j) + vvel(i+1,j))
                elseif (areafrac(i,j) > eps11) then
                   this_thck = thck_effective(i,j)
                   this_uvel = uvel(i,j)
                   this_vvel = vvel(i,j)
                else
                   this_thck = thck_effective(i+1,j)
                   this_uvel = uvel(i+1,j)
                   this_vvel = vvel(i+1,j)
                endif
                if (next_areafrac > eps11) then  ! take a weighted average between j and j+1
                   if (areafrac(i,j+1) > eps11 .and. areafrac(i+1,j+1) > eps11) then
                      next_thck = 0.5d0 * (thck_effective(i,j+1) + thck_effective(i+1,j+1))
                      next_uvel = 0.5d0 * (uvel(i,j+1) + uvel(i+1,j+1))
                      next_vvel = 0.5d0 * (vvel(i,j+1) + vvel(i+1,j+1))
                   elseif (areafrac(i,j+1) > eps11) then
                      next_thck = thck_effective(i,j+1)
                      next_uvel = uvel(i,j+1)
                      next_vvel = vvel(i,j+1)
                   else
                      next_thck = thck_effective(i+1,j+1)
                      next_uvel = uvel(i+1,j+1)
                      next_vvel = vvel(i+1,j+1)
                   endif
                   cf_thck(axis) = this_thck*wt_factor + next_thck*(1.0d0 - wt_factor)
                   cf_uvel(axis) = this_uvel*wt_factor + next_uvel*(1.0d0 - wt_factor)
                   cf_vvel(axis) = this_vvel*wt_factor + next_vvel*(1.0d0 - wt_factor)
                else   ! next_areafrac (at j+1) = 0; use values from this j
                   cf_thck(axis) = this_thck
                   cf_uvel(axis) = this_uvel
                   cf_vvel(axis) = this_vvel
                endif
             endif
          enddo
       endif
    enddo

    ! If this proc has a positive value of y, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_locy(axis), xout=cf_loc_ymax, xprocout=procnum)

    ! Broadcast the calvingMIP axis 2 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)


    axis = 4  ! index for Halbrane B (same as A except for positive x_intercept)
    x_intercept = 150.d3

    do i = nhalo, nx-nhalo
       if (abs(x0(i) - x_intercept) < eps11) then  ! E edge of cell lies on the vertical Halbrane profile
          cf_locx(axis) = x_intercept
          do j = nhalo, ny-nhalo
             this_areafrac = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
             next_areafrac = 0.5d0 * (areafrac(i,j+1) + areafrac(i+1,j+1))
             if (this_areafrac >= 0.5d0 .and. next_areafrac < 0.5d0) then
                ! CF lies between j and j+1
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locy(axis) = y1(j)*wt_factor + y1(j+1)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                ! The following logic allows for the possibility that one of the two neighbor cells is ice-free
                if (areafrac(i,j) > eps11 .and. areafrac(i+1,j) > eps11) then
                   this_thck = 0.5d0 * (thck_effective(i,j) + thck_effective(i+1,j))
                   this_uvel = 0.5d0 * (uvel(i,j) + uvel(i+1,j))
                   this_vvel = 0.5d0 * (vvel(i,j) + vvel(i+1,j))
                elseif (areafrac(i,j) > eps11) then
                   this_thck = thck_effective(i,j)
                   this_uvel = uvel(i,j)
                   this_vvel = vvel(i,j)
                else
                   this_thck = thck_effective(i+1,j)
                   this_uvel = uvel(i+1,j)
                   this_vvel = vvel(i+1,j)
                endif
                if (next_areafrac > eps11) then  ! take a weighted average between j and j+1
                   if (areafrac(i,j+1) > eps11 .and. areafrac(i+1,j+1) > eps11) then
                      next_thck = 0.5d0 * (thck_effective(i,j+1) + thck_effective(i+1,j+1))
                      next_uvel = 0.5d0 * (uvel(i,j+1) + uvel(i+1,j+1))
                      next_vvel = 0.5d0 * (vvel(i,j+1) + vvel(i+1,j+1))
                   elseif (areafrac(i,j+1) > eps11) then
                      next_thck = thck_effective(i,j+1)
                      next_uvel = uvel(i,j+1)
                      next_vvel = vvel(i,j+1)
                   else
                      next_thck = thck_effective(i+1,j+1)
                      next_uvel = uvel(i+1,j+1)
                      next_vvel = vvel(i+1,j+1)
                   endif
                   cf_thck(axis) = this_thck*wt_factor + next_thck*(1.0d0 - wt_factor)
                   cf_uvel(axis) = this_uvel*wt_factor + next_uvel*(1.0d0 - wt_factor)
                   cf_vvel(axis) = this_vvel*wt_factor + next_vvel*(1.0d0 - wt_factor)
                else   ! next_areafrac (at j+1) = 0; use values from this j
                   cf_thck(axis) = this_thck
                   cf_uvel(axis) = this_uvel
                   cf_vvel(axis) = this_vvel
                endif
             endif
          enddo
       endif
    enddo

    ! If this proc has a positive value of y, then broadcast the coordinates to all procs
    call parallel_reduce_maxloc(xin=cf_locy(axis), xout=cf_loc_ymax, xprocout=procnum)

    ! Broadcast the calvingMIP axis 4 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)

    
    axis = 6  ! index for Halbrane C (same as A except in the negative y direction)
    x_intercept = -150.d3

    do i = nhalo, nx-nhalo
       if (abs(x0(i) - x_intercept) < eps11) then  ! E edge of cell lies on the vertical Halbrane profile
          cf_locx(axis) = x_intercept
          do j = ny-nhalo, nhalo, -1
             this_areafrac = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
             next_areafrac = 0.5d0 * (areafrac(i,j-1) + areafrac(i+1,j-1))
             if (this_areafrac >= 0.5d0 .and. next_areafrac < 0.5d0) then
                ! CF lies between j and j-1
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locy(axis) = y1(j)*wt_factor + y1(j-1)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                ! The following logic allows for the possibility that one of the two neighbor cells is ice-free
                if (areafrac(i,j) > eps11 .and. areafrac(i+1,j) > eps11) then
                   this_thck = 0.5d0 * (thck_effective(i,j) + thck_effective(i+1,j))
                   this_uvel = 0.5d0 * (uvel(i,j) + uvel(i+1,j))
                   this_vvel = 0.5d0 * (vvel(i,j) + vvel(i+1,j))
                elseif (areafrac(i,j) > eps11) then
                   this_thck = thck_effective(i,j)
                   this_uvel = uvel(i,j)
                   this_vvel = vvel(i,j)
                else
                   this_thck = thck_effective(i+1,j)
                   this_uvel = uvel(i+1,j)
                   this_vvel = vvel(i+1,j)
                endif
                if (next_areafrac > eps11) then  ! take a weighted average between j and j-1
                   if (areafrac(i,j-1) > eps11 .and. areafrac(i+1,j-1) > eps11) then
                      next_thck = 0.5d0 * (thck_effective(i,j-1) + thck_effective(i+1,j-1))
                      next_uvel = 0.5d0 * (uvel(i,j-1) + uvel(i+1,j-1))
                      next_vvel = 0.5d0 * (vvel(i,j-1) + vvel(i+1,j-1))
                   elseif (areafrac(i,j-1) > eps11) then
                      next_thck = thck_effective(i,j-1)
                      next_uvel = uvel(i,j-1)
                      next_vvel = vvel(i,j-1)
                   else
                      next_thck = thck_effective(i+1,j-1)
                      next_uvel = uvel(i+1,j-1)
                      next_vvel = vvel(i+1,j-1)
                   endif
                   cf_thck(axis) = this_thck*wt_factor + next_thck*(1.0d0 - wt_factor)
                   cf_uvel(axis) = this_uvel*wt_factor + next_uvel*(1.0d0 - wt_factor)
                   cf_vvel(axis) = this_vvel*wt_factor + next_vvel*(1.0d0 - wt_factor)
                else   ! next_areafrac (at j-1) = 0; use values from this j
                   cf_thck(axis) = this_thck
                   cf_uvel(axis) = this_uvel
                   cf_vvel(axis) = this_vvel
                endif
             endif
          enddo
       endif
    enddo

    ! If this proc has a negative value of y, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_locy(axis), xout=cf_loc_ymin, xprocout=procnum)

    ! Broadcast the calvingMIP axis 6 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)

    
    axis = 8  ! index for Halbrane D (same as C except for positive x_intercept)
    x_intercept = 150.d3

    do i = nhalo, nx-nhalo
       if (abs(x0(i) - x_intercept) < eps11) then  ! E edge of cell lies on the vertical Halbrane profile
          cf_locx(axis) = x_intercept
          do j = ny-nhalo, nhalo, -1
             this_areafrac = 0.5d0 * (areafrac(i,j) + areafrac(i+1,j))
             next_areafrac = 0.5d0 * (areafrac(i,j-1) + areafrac(i+1,j-1))
             if (this_areafrac >= 0.5d0 .and. next_areafrac < 0.5d0) then
                ! CF lies between j and j-1
                wt_factor = (0.5d0 - next_areafrac) / (this_areafrac - next_areafrac)
                cf_locy(axis) = y1(j)*wt_factor + y1(j-1)*(1.0d0 - wt_factor)
                cf_radius(axis) = sqrt(cf_locx(axis)**2 + cf_locy(axis)**2)
                ! The following logic allows for the possibility that one of the two neighbor cells is ice-free
                if (areafrac(i,j) > eps11 .and. areafrac(i+1,j) > eps11) then
                   this_thck = 0.5d0 * (thck_effective(i,j) + thck_effective(i+1,j))
                   this_uvel = 0.5d0 * (uvel(i,j) + uvel(i+1,j))
                   this_vvel = 0.5d0 * (vvel(i,j) + vvel(i+1,j))
                elseif (areafrac(i,j) > eps11) then
                   this_thck = thck_effective(i,j)
                   this_uvel = uvel(i,j)
                   this_vvel = vvel(i,j)
                else
                   this_thck = thck_effective(i+1,j)
                   this_uvel = uvel(i+1,j)
                   this_vvel = vvel(i+1,j)
                endif
                if (next_areafrac > eps11) then  ! take a weighted average between j and j-1
                   if (areafrac(i,j-1) > eps11 .and. areafrac(i+1,j-1) > eps11) then
                      next_thck = 0.5d0 * (thck_effective(i,j-1) + thck_effective(i+1,j-1))
                      next_uvel = 0.5d0 * (uvel(i,j-1) + uvel(i+1,j-1))
                      next_vvel = 0.5d0 * (vvel(i,j-1) + vvel(i+1,j-1))
                   elseif (areafrac(i,j-1) > eps11) then
                      next_thck = thck_effective(i,j-1)
                      next_uvel = uvel(i,j-1)
                      next_vvel = vvel(i,j-1)
                   else
                      next_thck = thck_effective(i+1,j-1)
                      next_uvel = uvel(i+1,j-1)
                      next_vvel = vvel(i+1,j-1)
                   endif
                   cf_thck(axis) = this_thck*wt_factor + next_thck*(1.0d0 - wt_factor)
                   cf_uvel(axis) = this_uvel*wt_factor + next_uvel*(1.0d0 - wt_factor)
                   cf_vvel(axis) = this_vvel*wt_factor + next_vvel*(1.0d0 - wt_factor)
                else   ! next_areafrac (at j-1) = 0; use values from this j
                   cf_thck(axis) = this_thck
                   cf_uvel(axis) = this_uvel
                   cf_vvel(axis) = this_vvel
                endif
             endif
          enddo
       endif
    enddo

    ! If this proc has a negative value of y, then broadcast the coordinates to all procs
    call parallel_reduce_minloc(xin=cf_locy(axis), xout=cf_loc_ymin, xprocout=procnum)

    ! Broadcast the calvingMIP axis 8 output
    call broadcast(cf_locx(axis), proc=procnum)
    call broadcast(cf_locy(axis), proc=procnum)
    call broadcast(cf_radius(axis), proc=procnum)
    call broadcast(cf_thck(axis), proc=procnum)
    call broadcast(cf_uvel(axis), proc=procnum)
    call broadcast(cf_vvel(axis), proc=procnum)

    if (verbose_calvingmip .and. main_task) then
       write(iulog,*) ' '
       write(iulog,*) 'Thule domain: axis, CF location, radius (km), thck(m), uvel, vvel, speed (m/yr)'
       do axis = 1, 8
          speed = sqrt(cf_uvel(axis)**2 + cf_vvel(axis)**2)
          write(iulog,'(i4,7f15.8)') axis, cf_locx(axis)/1000.d0, cf_locy(axis)/1000.d0, &
               cf_radius(axis)/1000.d0, cf_thck(axis), cf_uvel(axis)*scyr, cf_vvel(axis)*scyr, speed*scyr
       enddo
    endif

  end subroutine locate_calving_front_thule

!---------------------------------------------------------------------------

  subroutine sum_over_quadrants(&
       nx,        ny,            &
       parallel,                 &
       field,     quadrant_sum)

    ! Integrate a field over each of 4 quadrants.
    ! This can be useful in idealized experiments like CalvingMIP to check for
    !  violations of reflectional or rotational symmetry.
    ! Note: These sums are not independent of processor count
    ! TODO: Make them reproducible, using quadrant masks?

    use cism_parallel, only: nhalo, parallel_global_sum_patch, parallel_globalindex, gather_var

    ! Input/output arguments

    integer, intent(in) :: &
         nx, ny                   ! number of local cells in x and y direction on input grid

    type(parallel_type), intent(in) :: parallel   ! info for parallel communication

    real(dp), dimension(nx,ny), intent(in) :: &
         field                    ! 2D input field

    real(dp), dimension(4), intent(out) :: &
         quadrant_sum             ! global sum over each of 4 quadrants

    logical, parameter :: check_asymmetry = .true.

    ! Local variables

    integer :: i, j
    integer :: ig, jg             ! i and j indices on the global grid
    integer :: nxg, nyg           ! dimensions of global domain
    integer :: nx2, ny2           ! nx/2 and ny/2 (if nx and ny are even)
                                  ! (nx-1)/2 and (ny-1)/2 (if nx and ny are odd)

    integer, dimension(nx,ny) :: &
         quadrant_mask            ! mask assigning each cell to a quadrent (1, 2, 3 or 4)

    real(dp), dimension(:,:), allocatable :: field_global

    real(dp) :: meanval, diff
    real(dp), parameter :: symmetry_tol = 1.0d-5    ! tolerance level for asymmetry


    ! Compute a mask that assigns each cell to one of 4 quadrants.
    ! Note: If nx or ny is odd, the middle row or column is excluded from the quadrant sums.

    nxg = parallel%global_ewn
    nyg = parallel%global_nsn

    if (mod(nxg,2) == 0) then   ! global_ewn is even
       nx2 = nxg/2
    else   ! nx is odd
       nx2 = (nxg-1)/2
    endif

    if (mod(nyg,2) == 0) then   ! global_nsn is even
       ny2 = nyg/2
    else   ! ny is odd
       ny2 = (nyg-1)/2
    endif

    quadrant_mask(:,:) = 0

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          call parallel_globalindex(i, j, ig, jg, parallel)
          if (ig > nx2) then
             if (jg > ny2) then   ! NE quadrant
                quadrant_mask(i,j) = 1
             else   ! jg <= ny2; SE quadrant
                quadrant_mask(i,j) = 4
             endif
          else   ! ig <= nx2
             if (jg > ny2) then   ! NW quadrant
                quadrant_mask(i,j) = 2
             else   ! jg <= ny2; SW quadrant
                quadrant_mask(i,j) = 3
             endif
          endif
       enddo
    enddo

    ! Compute the global sums
    ! Note: These sums are reproducible if reproducible_sums = .true.
    quadrant_sum(:) = parallel_global_sum_patch(field, 4, quadrant_mask, parallel)

    if (check_asymmetry) then

       call gather_var(field, field_global, parallel)

       if (main_task) then

          ! Identify asymmetries in reflection across the y-axis
          if (abs(quadrant_sum(1) + quadrant_sum(4) - quadrant_sum(2) - quadrant_sum(3)) > symmetry_tol) then
             do j = 1, nyg
                do i = 1, nx2
                   if (abs(field_global(i,j)) > eps11 .or. abs(field_global(nxg-i+1,j)) > eps11) then
                      meanval = 0.5d0 * (field_global(i,j) + field_global(nxg-i+1,j))
                      diff = abs(field_global(i,j) - field_global(nxg-i+1,j))
                      if (diff > meanval*symmetry_tol) then
                         write(iulog,*) 'Warning, y-reflection asymmetry: i, j, val(i,j), val(nxg-i+1,j), diff/mean:', &
                              i, j, field_global(i,j), field_global(nxg-i+1,j), diff/meanval
                      endif
                   endif
                enddo
             enddo
          endif

          ! Identify asymmetries in reflection across the x-axis
          if (abs(quadrant_sum(1) + quadrant_sum(2) - quadrant_sum(3) - quadrant_sum(4)) > symmetry_tol) then
             do j = 1, nyg
                do i = 1, nx2
                   if (abs(field_global(i,j)) > eps11 .or. abs(field_global(i,nyg-j+1)) > eps11) then
                      meanval = 0.5d0 * (field_global(i,j) + field_global(i,nyg-j+1))
                      diff = abs(field_global(i,j) - field_global(i,nyg-j+1))
                      if (diff > meanval*symmetry_tol) then
                         write(iulog,*) 'Warning, x-reflection asymmetry: i, j, val(i,j), val(i,nyg-j+1), diff/mean:', &
                              i, j, field_global(i,j), field_global(i,nyg-j+1), diff/meanval
                      endif
                   endif
                enddo
             enddo
          endif

          if (allocated(field_global)) deallocate(field_global)

       endif   ! main_task
    endif   ! check_asymmetry

  end subroutine sum_over_quadrants

!---------------------------------------------------------------------------

end module glissade_diagnostics

!---------------------------------------------------------------------------
