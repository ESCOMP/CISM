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
    use glissade_types
    use cism_parallel, only: this_rank, main_task, nhalo

    implicit none

    private
    public :: glissade_init_diag, glissade_write_diagnostics, &
         glissade_mass_balance_diagnostics, glissade_grounding_line_flux, &
         glissade_stress_tensor_eigenvalues, glissade_strain_rate_tensor_eigenvalues, &
         glissade_calvingmip_diag

    logical :: verbose_calvingmip = .false.

    logical, parameter :: verbose_diagnostics = .false.

  contains

!****************************************************************************

  subroutine glissade_init_diag (model)

    ! Initialize model diagnostics.
    ! (1) Set ndiag based on dt_diag.  (Diagnostics are written every ndiag steps.)
    ! (2) Find the local rank and indices of the global diagnostic point

    use glimmer_log, only: write_log, GM_DIAGNOSTIC
    use cism_parallel, only: parallel_localindex

    implicit none

    ! input/output arguments

    type(glissade_global_type), intent(inout) :: model    ! model instance

    ! local variables

    character(len=100) :: message

    !-----------------------------------------------------------------
    ! Given dt_diag, compute the interval ndiag of diagnostic output.
    ! (Output is written every ndiag timesteps.)
    ! NOTE: The ratio dt_diag/tinc is rounded to the nearest integer.
    !-----------------------------------------------------------------

    if (model%numerics%dt_diag > 0.0d0) then   ! dt_diag was specified in the config file; use it to compute ndiag

       ! Note: tinc and dt_diag have units of years, whereas dt has model timeunits
       model%numerics%ndiag = nint(model%numerics%dt_diag / model%numerics%tinc)
       model%numerics%ndiag = max(model%numerics%ndiag, 1)  ! cannot write more often than once per timestep

    endif

    !-----------------------------------------------------------------
    ! Find the local rank and indices of the global diagnostic point
    !-----------------------------------------------------------------

    call parallel_localindex(model%numerics%idiag,       model%numerics%jdiag,       &
                             model%numerics%idiag_local, model%numerics%jdiag_local, &
                             model%numerics%rdiag_local, model%parallel)

    if (main_task) then
       write(iulog,'(a26,2i6)') ' Global idiag, jdiag:     ',   &
                             model%numerics%idiag, model%numerics%jdiag
       write(iulog,'(a26,3i6)') ' Local idiag, jdiag, task:',   &
                             model%numerics%idiag_local,  &
                             model%numerics%jdiag_local,  &
                             model%numerics%rdiag_local
    endif

    if (main_task) then

       write(message,'(a25,2i6)') 'Global idiag, jdiag:     ',   &
                                   model%numerics%idiag, model%numerics%jdiag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,3i6)') 'Local idiag, jdiag, task:',   &
                                   model%numerics%idiag_local,  &
                                   model%numerics%jdiag_local,  &
                                   model%numerics%rdiag_local
       call write_log(trim(message), type = GM_DIAGNOSTIC)

    endif  ! main_task

    ! Broadcast from main task to all processors
    !TODO - Uncomment and make sure this does not cause problems
!    call broadcast(model%numerics%idiag_local)
!    call broadcast(model%numerics%jdiag_local)
!    call broadcast(model%numerics%rdiag_local)

  end subroutine glissade_init_diag

!****************************************************************************

  subroutine glissade_write_diagnostics (model,  time,    &
                                      tstep_count)

    ! Short driver subroutine to decide whether it's time to write diagnostics.
    ! If so, it calls glissade_write_diag.   

    ! input/output arguments

    ! Note: model is intent(inout) so that some global diagnostics can be computed below
    type(glissade_global_type), intent(inout) :: model    ! model instance

    real(dp), intent(in) :: time          ! current time in years

    integer,  intent(in) :: tstep_count   ! current timestep

    ! local arguments

    ! debug
    if (main_task .and. verbose_diagnostics) then
       write(iulog,*) '	'
       write(iulog,*) 'In glissade_write_diagnostics'
       write(iulog,*) 'time =', time
       write(iulog,*) 'dt_diag =', model%numerics%dt_diag
       write(iulog,*) 'ndiag =', model%numerics%ndiag
       write(iulog,*) 'tstep_count =', tstep_count
    endif

    if (model%numerics%ndiag > 0) then

       if (mod(tstep_count, model%numerics%ndiag) == 0)  then    ! time to write
          call glissade_write_diag(model, time)
       endif

    endif    ! ndiag > 0

  end subroutine glissade_write_diagnostics

!****************************************************************************

  subroutine glissade_write_diag (model,       time)

    ! Write global diagnostics
    ! Also write local diagnostics for a selected grid cell
 
    use glimmer_paramets, only: unphys_val
    use glimmer_physcon, only: scyr, rhoi, shci
    use glimmer_log
    use cism_parallel, only: lhalo, uhalo, parallel_type, broadcast, &
         parallel_globalindex, &
         parallel_global_sum, parallel_global_sum_patch, &
         parallel_reduce_max, parallel_reduce_maxloc, parallel_reduce_minloc, &
         parallel_is_zero
    use glissade_utils, only: glissade_usrf_to_thck, glissade_rms_error

    implicit none
 
    ! input/output arguments

    ! Note: model is intent(inout) so that some global diagnostics can be computed here 
    type(glissade_global_type), intent(inout) :: model ! model instance

    real(dp),  intent(in) :: time                   ! current time in years

    ! local variables

    real(dp) ::                         &
         minthck,                       &    ! ice thickness threshold (m) for global diagnostics
         tot_area,                      &    ! total ice area (m^2)
         tot_area_ground,               &    ! total area of grounded ice (m^2)
         tot_area_float,                &    ! total area of floating ice (m^2)
         area_cell,                     &    ! cell area
         tot_volume,                    &    ! total ice volume (m^3)
         tot_volume_above_flotation,    &    ! total ice volume above flotation (m^3)
         tot_mass,                      &    ! total ice mass (kg)
         tot_mass_above_flotation,      &    ! total ice mass above flotation (kg)
         tot_area_ice_caps,             &    ! total area of disconnected ice caps (m^2)
         tot_vol_ice_caps,              &    ! total volume of disconnected ice caps (m^3)
         tot_bmlt_float,                &    ! total basal melt rate for floating ice (kg/s)
         tot_bmlt_float_target,         &    ! total target basal melt rate for floating ice (kg/s)
         thck_floating,                 &    ! thickness of floating ice
         thck_above_flotation,          &    ! thickness above flotation
         tot_energy,                    &    ! total ice energy (J)
         tot_smb_flux,                  &    ! total surface mass balance flux (kg/s)
         tot_bmb_flux,                  &    ! total basal mass balance flux (kg/s)
         tot_calving_flux,              &    ! total calving flux (kg/s)
         tot_latmelt_flux,              &    ! total latmelt flux (kg/s)
         tot_removal_flux,              &    ! total removal flux (kg/s)
         tot_gl_flux,                   &    ! total grounding line flux (kg/s)
         tot_acab,                      &    ! total surface accumulation/ablation rate (m^3/yr)
         tot_bmlt,                      &    ! total basal melt rate (m^3/yr)
         tot_calving,                   &    ! total calving rate (m^3/yr)
         tot_latmelt,                   &    ! total lateral melt rate (m^3/yr)
         tot_removal,                   &    ! total removal rate (m^3/yr)
         tot_dmass_dt,                  &    ! rate of change of total mass (kg/s)
         err_dmass_dt,                  &    ! mass conservation error (kg/s)
                                             ! given by dmass_dt - (tot_acab - tot_bmlt - tot_calving)
         mean_thck,                     &    ! mean ice thickness (m)
         mean_temp,                     &    ! mean ice temperature (deg C)
         mean_acab,                     &    ! mean surface accumulation/ablation rate (m/yr)
         mean_bmlt,                     &    ! mean basal melt (m/yr)
         mean_calving,                  &    ! mean calving (m/yr)
         mean_latmelt,                  &    ! mean latmelt (m/yr)
         mean_removal,                  &    ! mean removal (m/yr)
         max_thck, max_thck_global,     &    ! max ice thickness (m)
         max_temp, max_temp_global,     &    ! max ice temperature (deg C)
         min_temp, min_temp_global,     &    ! min ice temperature (deg C)
         max_bmlt,                      &    ! max basal melt rate (m/yr)
         max_bmlt_global,               &
         max_bmlt_ground,               &    ! max basal melt rate, grounded ice (m/yr)
         max_bmlt_ground_global,        &
         max_spd_sfc, max_spd_sfc_global, &  ! max surface ice speed (m/yr)
         max_spd_bas, max_spd_bas_global, &  ! max basal ice speed (m/yr)
         spd,                           &    ! speed
         thck_diag, usrf_diag,          &    ! local column diagnostics
         topg_diag, relx_diag,          &    
         load_diag,                     &
         artm_diag, acab_diag,          &
         bmlt_diag, bwat_diag,          &
         bheatflx_diag,                 &
         top_age_diag, bot_age_diag,    &
         level,                         &
         rmse_thck, rmse_velo,          &    ! rms errors for ice thickness (m) and surface speed (m/yr)
         factor                              ! unit conversion factor

    integer, dimension(model%general%ewn,model%general%nsn) ::  &
         ice_mask,                 & ! = 1 where ice is present with thck > minthck, else = 0
         floating_mask,            & ! = 1 where ice is present and floating, else = 0
         grounded_mask,            & ! = 1 where ice is present and grounded, else = 0
         ice_cap_mask,             & ! = 1 where an ice cap is present, else = 0
         glacier_ice_mask            ! = 1 where glacier ice is present, initially and/or currently

    integer, dimension(model%general%ewn-1,model%general%nsn-1) ::  &
         stag_ice_mask               ! staggered mask; = 1 if ice_mask = 1 for any of the four neighbors

    real(dp), dimension(model%general%upn) ::  &
         temp_diag,                     &    ! Note: sfc temp not included if temps are staggered
                                             !       (use artm instead)
         spd_diag                            ! speed (m/yr)

    real(dp), dimension(model%general%upn-1) ::  &
         age_diag                            ! ice age (yr)

    real(dp), dimension(model%lithot%nlayer) ::  &
         lithtemp_diag                       ! lithosphere column diagnostics

    ! glacier diagnostics
    real(dp) :: &
         tot_glc_area_init, tot_glc_area,     & ! total glacier area, initial and current (km^2)
         tot_glc_volume_init, tot_glc_volume, & ! total glacier volume, initial and current (km^3)
         tot_glc_area_init_extent,            & ! glacier area summed over the initial extent (km^2)
         tot_glc_volume_init_extent,          & ! glacier volume summed over the initial extent (km^3)
         tot_glc_area_target,                 & ! target glacier area for inversion (km^2)
         tot_glc_volume_target,               & ! target glacier volume for inversion (km^3)
         glc_rmse_thck,                       & ! root mean square value of thck - thck_target
         glc_rmse_thck_init_extent              ! as above, but within initial extent

    integer :: &
         count_area, count_volume               ! number of glaciers with nonzero area and volume

    integer :: &
         i, j, k, nb, ng,                   &
         ktop, kbed,                        &
         imax, imin,                        &
         jmax, jmin,                        &
         kmax, kmin,                        &
         imax_global, imin_global,          &
         jmax_global, jmin_global,          &
         kmax_global, kmin_global,          &
         procnum,                           &
         ewn, nsn, upn,                     &    ! model%numerics%ewn, etc.
         nlith,                             &    ! model%lithot%nlayer
         velo_ew_ubound, velo_ns_ubound          ! upper bounds for velocity variables

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         volume_above_flotation,& ! ice volue above flotation (m^3)
         thck_obs                 ! observed ice thickness (m), derived from usrf_obs and topg

    real(dp), dimension(model%general%ewn-1, model%general%nsn-1) ::  &
         velo_sfc               ! surface ice speed (m/s)

    real(dp), dimension(model%general%upn, model%general%ewn, model%general%nsn) ::  &
         local_energy           ! internal energy (J) per layer in a column

    character(len=100) :: message
    
    ! Note: cell_area is copied here from model%geometry%cell_area
    ! cell_area = dew*dns by default; optionally scaled to account for grid distortion
    real(dp), dimension(model%general%ewn,model%general%nsn) :: &
         cell_area     ! grid cell areas (scaled model units); diagnostic only

    real(dp), parameter ::   &
         eps = eps11             ! small threshold for diagnostics

    type(parallel_type) :: parallel       ! info for parallel communication

    parallel = model%parallel
    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn

    ! Set cell_area = model%geometry%cell_area
    ! Note: By default, cell_area = dew*dns
    !       For diagnostics, however, we may want to correct for grid distortions,
    !        giving a better estimate of the true ice area and volume.
    !       In this case, model%geometry%cell_area is corrected at initialization.
    !       It is used only for diagnostics. In the dynamics, each cell is a rectangle of area dew*dns.
    !       Using the corrected value here will give a conservation error (total dmass_dt > 0)
    !        in the diagnostics, because horizontal transport does not account for area factors.
    !        Horizontal transport conserves mass only under the assumption of rectangular grid cells.
    cell_area = model%geometry%cell_area

    nlith = model%lithot%nlayer

    if (uhalo > 0) then
       velo_ns_ubound = nsn-uhalo
       velo_ew_ubound = ewn-uhalo
    else
       ! for uhalo==0 (as is the case for the glide dycore), the velocity grid has one less
       ! point than the main grid, so we need to subtract one to avoid out-of-bounds problems
       velo_ns_ubound = nsn-uhalo-1
       velo_ew_ubound = ewn-uhalo-1
    end if

    ! Set the minimum ice thickness for including cells in diagnostics
    if (model%options%diag_minthck == DIAG_MINTHCK_ZERO) then
       minthck = eps  ! slightly > 0
    elseif (model%options%diag_minthck == DIAG_MINTHCK_THKLIM) then
       ! Note: If the user specifies model%numerics%thklim = 1.0 m in the config file,
       !       then thklim is reset at initialization to a number slightly less than 1.0.
       !       This protects the diagnostics below against rounding errors for cells
       !       with thck very close to minthck.
       minthck = model%numerics%thklim
    endif

    !-----------------------------------------------------------------
    ! Compute some masks that are useful for diagnostics
    !-----------------------------------------------------------------

    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > minthck) then
             ice_mask(i,j) = 1
             if (model%geometry%topg(i,j) - model%climate%eus < (-rhoi/rhoo)*model%geometry%thck(i,j)) then
                floating_mask(i,j) = 1
                grounded_mask(i,j) = 0
             else
                floating_mask(i,j) = 0
                grounded_mask(i,j) = 1
             endif
          else
             ice_mask(i,j) = 0
             floating_mask(i,j) = 0
             grounded_mask(i,j) = 0
          endif
       enddo
    enddo

    do j = 1, nsn-1
       do i = 1, ewn-1
          if (ice_mask(i,j+1) == 1 .or. ice_mask(i+1,j+1) == 1 .or.  &
              ice_mask(i,j)   == 1 .or. ice_mask(i+1,j)   == 1) then
             stag_ice_mask(i,j) = 1
          else
             stag_ice_mask(i,j) = 0
          endif
       enddo
    enddo

    ice_cap_mask = model%geometry%ice_cap_mask

    !-----------------------------------------------------------------
    ! Compute and write global diagnostics
    !-----------------------------------------------------------------
 
    call write_log('----------------------------------------------------------')
    call write_log(' ')
    write(message,'(a25,f24.16)') 'Diagnostic output, time =', time
    call write_log(trim(message), type = GM_DIAGNOSTIC)
    call write_log(' ')

    ! total ice area (m^2)
    tot_area = parallel_global_sum(cell_area, parallel, ice_mask)
    tot_area_float = parallel_global_sum(cell_area, parallel, floating_mask)
    tot_area_ground = parallel_global_sum(cell_area, parallel, grounded_mask)

    ! total ice volume (m^3)
    tot_volume = parallel_global_sum(model%geometry%thck*cell_area, parallel, ice_mask)

    ! total ice mass (kg)
    tot_mass = tot_volume * rhoi

    ! total ice volume above flotation (m^3)
    volume_above_flotation = 0.0d0
    do j = 1, nsn
       do i = 1, ewn
          if (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) then
             if (model%geometry%topg(i,j) - model%climate%eus < 0.0d0) then  ! grounded below sea level
                thck_floating = (-rhoo/rhoi) * (model%geometry%topg(i,j) - model%climate%eus)  ! exactly floating
                volume_above_flotation(i,j) = (model%geometry%thck(i,j) - thck_floating) * cell_area(i,j)
             else   ! grounded above sea level
                volume_above_flotation(i,j) = model%geometry%thck(i,j) * cell_area(i,j)
             endif
          endif
       enddo
    enddo

    tot_volume_above_flotation = parallel_global_sum(volume_above_flotation, parallel)

    ! total ice mass above flotation (kg)
    tot_mass_above_flotation = tot_volume_above_flotation * rhoi

    ! ice cap area and volume
    tot_area_ice_caps = parallel_global_sum(cell_area, parallel, ice_cap_mask)
    tot_vol_ice_caps = parallel_global_sum(model%geometry%thck*cell_area, parallel, ice_cap_mask)

    ! basal melting for floating ice (kg/s)
    tot_bmlt_float = &
         parallel_global_sum(model%basal_melt%bmlt_float*rhoi*cell_area, parallel)
    tot_bmlt_float_target = &
         parallel_global_sum(model%basal_melt%bmlt_float_target*rhoi*cell_area, parallel)

    ! total ice energy relative to T = 0 deg C (J)
    local_energy = 0.0d0
    if (size(model%temper%temp,1) == upn+1) then  ! temps are staggered in vertical, located at layer centers
       ! in this case, the data in the upn layer remains = 0
       do j = 1, nsn
          do i = 1, ewn
             if (ice_mask(i,j) == 1) then
                do k = 1, upn-1   ! (upn-1) layers
                   local_energy(k,i,j) = &
                        model%geometry%thck(i,j) * model%temper%temp(k,i,j) * cell_area(i,j)   &
                        * (model%numerics%sigma(k+1) - model%numerics%sigma(k))
                enddo
             endif
          enddo
       enddo
    else   ! temps are unstaggered in vertical, located at layer interfaces
       do j = 1, nsn
          do i = 1, ewn
             if (ice_mask(i,j) == 1) then
                k = 1   ! top interface; assign this temperature to the top half layer
                local_energy(k,i,j) = &
                     model%geometry%thck(i,j) * model%temper%temp(k,i,j) * cell_area(i,j)    &
                     * 0.5d0 * model%numerics%sigma(k+1)
                do k = 2, upn-1   ! interior layers
                   local_energy(k,i,j) = &
                        model%geometry%thck(i,j) * model%temper%temp(k,i,j) * cell_area(i,j)  &
                        * 0.5d0*(model%numerics%sigma(k+1) - model%numerics%sigma(k-1))
                enddo
                k = upn   ! bottom interface; assign this temperature to the bottom half layer
                local_energy(k,i,j) = &
                     model%geometry%thck(i,j) * model%temper%temp(k,i,j) * cell_area(i,j)  &
                     * 0.5d0 * (1.0d0 - model%numerics%sigma(k-1))
             endif
          enddo
       enddo
    endif
    local_energy = local_energy * rhoi * shci
    tot_energy = parallel_global_sum(local_energy, parallel, ice_mask)

    ! mean thickness
    if (tot_area > eps) then
       mean_thck = tot_volume/tot_area
    else
       mean_thck = 0.d0
    endif

    ! mean temperature
    if (tot_volume > eps) then
       mean_temp = tot_energy/ (rhoi*shci*tot_volume)
    else
       mean_temp = 0.d0
    endif
 
    ! copy some global scalars to the scalars derived type
    ! Note: These have SI units (e.g, m^2 for area, m^3 for volume, kg/s for melt rates)

    model%scalars%iarea  = tot_area
    model%scalars%iareag = tot_area_ground
    model%scalars%iareaf = tot_area_float
    model%scalars%ivol   = tot_volume
    model%scalars%ivol_above_flotation = tot_volume_above_flotation
    model%scalars%imass  = tot_mass
    model%scalars%imass_above_flotation = tot_mass_above_flotation
    model%scalars%icap_area  = tot_area_ice_caps
    model%scalars%icap_vol  = tot_vol_ice_caps
    model%scalars%total_bmlt_float = tot_bmlt_float
    model%scalars%total_bmlt_float_target = tot_bmlt_float_target

    ! Optionally, compute some basin-scale scalars, also written to the scalars derived type

    if (model%ocean_data%nbasin > 1) then
       model%scalars%iarea_basin(:)  = &
            parallel_global_sum_patch(cell_area*ice_mask, model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)
       model%scalars%iareag_basin(:) = &
            parallel_global_sum_patch(cell_area*grounded_mask, model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)
       model%scalars%iareaf_basin(:) = &
            parallel_global_sum_patch(cell_area*floating_mask, model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)
       model%scalars%ivol_basin(:)   = &
            parallel_global_sum_patch(cell_area*ice_mask*model%geometry%thck, model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)
       model%scalars%ivol_above_flotation_basin(:) = &
            parallel_global_sum_patch(volume_above_flotation, model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)
       model%scalars%imass_basin(:) = model%scalars%ivol_basin(:)*rhoi
       model%scalars%imass_above_flotation_basin(:) = model%scalars%ivol_above_flotation_basin(:)*rhoi
       model%scalars%icap_area_basin(:)  = &
            parallel_global_sum_patch(cell_area*ice_cap_mask, model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)
       model%scalars%icap_vol_basin(:)  = &
            parallel_global_sum_patch(cell_area*ice_cap_mask*model%geometry%thck, &
            model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)
       model%scalars%bmlt_float_basin(:)  = &
            parallel_global_sum_patch(model%basal_melt%bmlt_float*rhoi*cell_area, &
            model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)
       model%scalars%bmlt_float_target_basin(:)  = &
            parallel_global_sum_patch(model%basal_melt%bmlt_float_target*rhoi*cell_area, &
            model%ocean_data%nbasin, model%ocean_data%basin_number, parallel)

       ! Optionally, write output to a specific basin with an applied thermal forcing anomaly
       if (main_task) then
          nb = model%ocean_data%thermal_forcing_anomaly_basin
          if (nb >= 1 .and. nb <= model%ocean_data%nbasin) then
             write(iulog,*) 'Diagnostics for basin', nb
             write(iulog,*) 'iarea, iareag, iareaf (km^2):', &
                  model%scalars%iarea_basin(nb)/1.0d6, model%scalars%iareag_basin(nb)/1.0d6, model%scalars%iareaf_basin(nb)/1.0d6
             write(iulog,*) 'ivol, ivol_above_flotation (km^3):', &
                  model%scalars%ivol_basin(nb)/1.0d9, model%scalars%ivol_above_flotation_basin(nb)/1.0d9
             write(iulog,*) 'imass, imass_above_flotation (Gt):', &
                  model%scalars%imass_basin(nb)/1.0d12, model%scalars%imass_above_flotation_basin(nb)/1.0d12
             if (.not.model%options%remove_ice_caps)  then
                write(iulog,*) '   ice cap area (km^2):', model%scalars%icap_area_basin(nb)/1.0d6
                write(iulog,*) '   ice cap vol  (km^3):', model%scalars%icap_vol_basin(nb)/1.0d9
             endif
          endif
       endif

    endif   ! nbasin > 1

    ! Compute a global mass budget and check mass conservation

    ! total surface accumulation/ablation rate (m^3/yr ice)
    tot_acab = parallel_global_sum(model%climate%acab_applied*cell_area, parallel)
    tot_acab = tot_acab * scyr                ! convert from m^3/s to m^3/yr

    ! total surface mass balance flux (kg/s)
    tot_smb_flux = tot_acab * rhoi / scyr     ! convert m^3/yr to kg/s

    ! mean accumulation/ablation rate (m/yr)
    ! Note: This will be only approximate if some ice has melted completely during the time step
    if (tot_area > eps) then
       mean_acab = tot_acab/tot_area    ! divide by total area to get m/yr
    else
       mean_acab = 0.d0
    endif

    ! total basal melting rate (positive for ice loss)
    tot_bmlt = parallel_global_sum(model%basal_melt%bmlt_applied*cell_area, parallel)
    tot_bmlt = tot_bmlt * scyr                ! convert from m^3/s to m^3/yr

    ! total basal mass balance (kg/s, positive for freeze-on, negative for melt)
    tot_bmb_flux = -tot_bmlt * rhoi / scyr   ! convert m^3/yr to kg/s

    ! mean basal melt rate (m/yr)
    ! Note: This will be only approximate if some ice has melted completely during the time step
    if (tot_area > eps) then
       mean_bmlt = tot_bmlt/tot_area    ! divide by total area to get m/yr
    else
       mean_bmlt = 0.d0
    endif

    ! total calving rate (m^3/yr ice)
    ! Note: calving%calving_rate has units of m/yr ice
    tot_calving = parallel_global_sum(model%calving%calving_rate*cell_area, parallel)

    ! total calving mass balance flux (kg/s, negative for ice loss by calving)
    tot_calving_flux = -tot_calving * rhoi / scyr   ! convert m^3/yr to kg/s

    ! mean calving rate (m/yr)
    ! Note: This will be only approximate if some ice has calved completely during the time step
    if (tot_area > eps) then
       mean_calving = tot_calving/tot_area    ! divide by total area to get m/yr
    else
       mean_calving = 0.d0
    endif

    ! total lateral melt rate (m^3/yr ice)
    ! Note: latmelt_rate has units of m/yr ice
    tot_latmelt = parallel_global_sum(model%lateral_melt%latmelt_rate*cell_area, parallel)

    ! total lateral melt mass balance flux (kg/s, negative for ice loss by melting)
    tot_latmelt_flux = -tot_latmelt * rhoi/scyr    ! convert m^3/yr to kg/s

    ! mean lateral melt rate (m/yr)
    ! Note: This will be only approximate if some ice has melted completely during the time step
    if (tot_area > eps) then
       mean_latmelt = tot_latmelt/tot_area    ! divide by total area to get m/yr
    else
       mean_latmelt = 0.d0
    endif

    ! total ice removal rate (m^3/yr ice)
    ! Note: geometry%removal_rate has units of m/yr ice
    tot_removal = parallel_global_sum(model%geometry%removal_rate*cell_area, parallel)

    ! total removal mass balance flux (kg/s, negative for ice removed)
    tot_removal_flux = -tot_removal * rhoi / scyr   ! convert m^3/yr to kg/s

    ! mean removal rate (m/yr)
    if (tot_area > eps) then
       mean_removal = tot_removal/tot_area    ! divide by total area to get m/yr
    else
       mean_removal = 0.d0
    endif

    ! total grounding line mass balance flux (< 0 by definition)
    ! Note: At this point, gl_flux_east and gl_flux_north are already dimensionalized in kg/m/s,
    !       so tot_gl_flux will have units of kg/s
    tot_gl_flux = parallel_global_sum(abs(model%mass_flux%gl_flux_east)  * model%numerics%dns  &
                                    + abs(model%mass_flux%gl_flux_north) * model%numerics%dew, &
                                    parallel)
    tot_gl_flux = (-1.0d0)*tot_gl_flux   ! negative by definition

    ! total rate of change of ice mass (kg/s)
    ! Note: dthck_dt has units of m/s
    !TODO: Change this calculation to use the total mass of ice in the global domain
    !       in successive time steps, instead of summing over dthck_dt.
    !      Note that dthck_dt does not account for global outflow fluxes (i.e., removal of ice
    !       near the global boundary in halo updates).
    tot_dmass_dt = parallel_global_sum(model%geometry%dthck_dt*cell_area, parallel)
    tot_dmass_dt = tot_dmass_dt * rhoi              ! convert to kg/s

    ! mass conservation error
    ! Note: For most runs, this should be close to zero.

    err_dmass_dt = tot_dmass_dt - &
         (tot_smb_flux + tot_bmb_flux + tot_calving_flux + tot_latmelt_flux + tot_removal_flux)

    ! copy some global scalars to the mass_flux derived type
    ! Note: These have SI units (e.g, m^2 for area, m^3 for volume)
    model%mass_flux%total_smb_flux = tot_smb_flux
    model%mass_flux%total_bmb_flux = tot_bmb_flux
    model%mass_flux%total_calving_flux = tot_calving_flux
    model%mass_flux%total_latmelt_flux = tot_latmelt_flux
    model%mass_flux%total_removal_flux = tot_removal_flux

    ! Note: The total mass budget consists of the five terms above.
    !       The GL flux is not part of the mass budget, since any ice fluxed across the GL
    !        is not added or lost, but simply changes from grounded to floating.
    model%mass_flux%total_gl_flux = tot_gl_flux

    ! write global sums and means to log file

    write(message,'(a25,e24.16)') 'Total ice area (km^2)    ',   &
                                   tot_area*1.0d-6           ! convert to km^2
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Grounded ice area (km^2) ',   &
                                   tot_area_ground*1.0d-6    ! convert to km^2
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Floating ice area (km^2) ',   &
                                   tot_area_float*1.0d-6     ! convert to km^2
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Total ice volume (km^3)  ',   &
                                   tot_volume*1.0d-9         ! convert to km^3
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    if (model%options%dm_dt_diag == DM_DT_DIAG_KG_S) then

       write(message,'(a25,e24.16)') 'Total ice mass (kg)      ',   &
                                      tot_mass                  ! kg
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'Mass above flotation (kg)',   &
                                      tot_mass_above_flotation  ! kg
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'Total ice energy (J)     ',   &
                                      tot_energy                ! J
       call write_log(trim(message), type = GM_DIAGNOSTIC)

    elseif (model%options%dm_dt_diag == DM_DT_DIAG_GT_Y) then

       write(message,'(a25,e24.16)') 'Total ice mass (Gt)      ',   &
                                      tot_mass * 1.0d-12        ! Gt
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'Mass above flotation (Gt)',   &
                      tot_mass_above_flotation * 1.0d-12        ! Gt
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'Total ice energy (GJ)     ',  &
                                    tot_energy * 1.0d-9         ! GJ
       call write_log(trim(message), type = GM_DIAGNOSTIC)

    endif  ! dm_dt_diag

    if (.not.model%options%remove_ice_caps .and. .not.model%options%enable_glaciers) then
       write(message,'(a25,e24.16)') 'Ice cap area (km^2)       ',   &
                                      tot_area_ice_caps*1.0d-6  ! convert to km^2
       call write_log(trim(message), type = GM_DIAGNOSTIC)
       write(message,'(a25,e24.16)') 'Ice cap volume (km^3)     ',   &
                                      tot_vol_ice_caps*1.0d-9   ! convert to km^3
       call write_log(trim(message), type = GM_DIAGNOSTIC)
    endif

    if (model%options%dm_dt_diag == DM_DT_DIAG_KG_S) then

       write(message,'(a25,e24.16)') 'Total SMB flux (kg/s)    ', tot_smb_flux
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'Total BMB flux (kg/s)    ', tot_bmb_flux
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'Total calving flux (kg/s)', tot_calving_flux
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       if (abs(tot_latmelt_flux) > eps11) then
          write(message,'(a25,e24.16)') 'Total latmelt flux (kg/s)', tot_latmelt_flux
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       endif

       if (abs(tot_removal_flux) > eps11) then
          write(message,'(a25,e24.16)') 'Total removal flux (kg/s)', tot_removal_flux
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       endif

       if (tot_bmlt_float > eps11) then
          write(message,'(a25,e24.16)') 'Total bmlt_float (kg/s)  ', tot_bmlt_float
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       endif

       if (tot_bmlt_float_target > eps11) then
          write(message,'(a25,e24.16)') 'Target bmlt_float (kg/s)', tot_bmlt_float_target
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       endif

       write(message,'(a25,e24.16)') 'Total dmass/dt (kg/s)    ', tot_dmass_dt
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'dmass/dt error (kg/s)    ', err_dmass_dt
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'Total gr line flux (kg/s)', tot_gl_flux
       call write_log(trim(message), type = GM_DIAGNOSTIC)

    elseif (model%options%dm_dt_diag == DM_DT_DIAG_GT_Y) then

       factor = scyr / 1.0d12

       write(message,'(a25,e24.16)') 'Total SMB flux (Gt/y)    ', tot_smb_flux * factor
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'Total BMB flux (Gt/y)    ', tot_bmb_flux * factor
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'Total calving flux (Gt/y)', tot_calving_flux * factor
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       if (abs(tot_latmelt_flux) > eps11) then
          write(message,'(a25,e24.16)') 'Total latmelt flux (Gt/y)', tot_latmelt_flux * factor
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       endif

       if (abs(tot_removal_flux) > eps11) then
          write(message,'(a25,e24.16)') 'Total removal flux (Gt/y)', tot_removal_flux * factor
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       endif

       if (tot_bmlt_float > eps11) then
          write(message,'(a25,e24.16)') 'Total bmlt_float (Gt/y)  ', tot_bmlt_float * factor
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       endif

       if (tot_bmlt_float_target > eps11) then
          write(message,'(a25,e24.16)') 'Target bmlt_float (Gt/y) ', tot_bmlt_float_target * factor
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       endif

       write(message,'(a25,e24.16)') 'Total dmass/dt (Gt/y)    ', tot_dmass_dt * factor
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'dmass/dt error (Gt/y)    ', err_dmass_dt * factor
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,e24.16)') 'Total gr line flux (Gt/y)', tot_gl_flux * factor
       call write_log(trim(message), type = GM_DIAGNOSTIC)

    endif   ! dm_dt_diag

!       write(message,'(a25,e24.16)') 'Mean accum/ablat (m/yr)  ', mean_acab
!       call write_log(trim(message), type = GM_DIAGNOSTIC)

!       write(message,'(a25,e24.16)') 'Mean basal melt (m/yr)   ', mean_bmlt
!       call write_log(trim(message), type = GM_DIAGNOSTIC)

!       write(message,'(a25,e24.16)') 'Mean calving (m/yr)      ', mean_calving
!       call write_log(trim(message), type = GM_DIAGNOSTIC)

!       write(message,'(a25,e24.16)') 'Mean latmelt (m/yr)      ', mean_latmelt
!       call write_log(trim(message), type = GM_DIAGNOSTIC)

!       write(message,'(a25,e24.16)') 'Mean removal (m/yr)      ', mean_removal
!       call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,f24.16)') 'Mean thickness (m)       ', mean_thck
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,f24.16)') 'Mean temperature (C)     ', mean_temp
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! Find various global maxes and mins

    ! max thickness

    imax = 0
    jmax = 0
    max_thck = unphys_val   ! = an arbitrary large negative number
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%geometry%thck(i,j) > max_thck) then
             max_thck = model%geometry%thck(i,j)
             imax = i
             jmax = j
          endif
       enddo
    enddo

    imax_global = 0
    jmax_global = 0
    max_thck_global = parallel_reduce_max(max_thck)
    if (max_thck == max_thck_global) then  ! max_thck lives on this processor
       imax_global = (imax - lhalo) + parallel%global_col_offset
       jmax_global = (jmax - lhalo) + parallel%global_row_offset
    endif
    imax_global = parallel_reduce_max(imax_global)
    jmax_global = parallel_reduce_max(jmax_global)

    write(message,'(a25,f24.16,2i6)') 'Max thickness (m), i, j  ',   &
                                       max_thck_global, imax_global, jmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max temperature

    ktop = lbound(model%temper%temp,1)
    kbed = ubound(model%temper%temp,1)

    imax = 0
    jmax = 0
    kmax = 0
    max_temp =  unphys_val
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (ice_mask(i,j) == 1) then
            do k = ktop, kbed
                if (model%temper%temp(k,i,j) > max_temp) then
                   max_temp = model%temper%temp(k,i,j)
                   imax = i
                   jmax = j
                   kmax = k
                endif
             enddo
          endif
       enddo
    enddo

    call parallel_reduce_maxloc(xin=max_temp, xout=max_temp_global, xprocout=procnum)
    call parallel_globalindex(imax, jmax, imax_global, jmax_global, parallel)
    kmax_global = kmax
    call broadcast(imax_global, proc = procnum)
    call broadcast(jmax_global, proc = procnum)
    call broadcast(kmax_global, proc = procnum)

    write(message,'(a25,f24.16,3i6)') 'Max temperature, i, j, k ',   &
                    max_temp_global, imax_global, jmax_global, kmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! min temperature

    imin = 0
    jmin = 0
    kmin = 0
    min_temp = 999.d0  ! arbitrary large positive number
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (ice_mask(i,j) == 1) then
             do k = ktop, kbed
                if (model%temper%temp(k,i,j) < min_temp) then
                   min_temp = model%temper%temp(k,i,j)
                   imin = i
                   jmin = j
                   kmin = k
                endif
             enddo
          endif
       enddo
    enddo

    call parallel_reduce_minloc(xin=min_temp, xout=min_temp_global, xprocout=procnum)
    call parallel_globalindex(imin, jmin, imin_global, jmin_global, parallel)
    kmin_global = kmin
    call broadcast(imin_global, proc = procnum)
    call broadcast(jmin_global, proc = procnum)
    call broadcast(kmin_global, proc = procnum)

    write(message,'(a25,f24.16,3i6)') 'Min temperature, i, j, k ',   &
                    min_temp_global, imin_global, jmin_global, kmin_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max applied basal melt rate
    ! Usually, this will be for floating ice, if floating ice is present
    imax = 0
    jmax = 0
    max_bmlt = unphys_val

    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%basal_melt%bmlt_applied(i,j) > max_bmlt) then
             max_bmlt = model%basal_melt%bmlt_applied(i,j)
             imax = i
             jmax = j
          endif
       enddo
    enddo

    call parallel_reduce_maxloc(xin=max_bmlt, xout=max_bmlt_global, xprocout=procnum)
    call parallel_globalindex(imax, jmax, imax_global, jmax_global, parallel)
    call broadcast(imax_global, proc = procnum)
    call broadcast(jmax_global, proc = procnum)

    ! Write to diagnostics only if nonzero

    if (abs(max_bmlt_global*scyr) > eps) then
       write(message,'(a25,f24.16,2i6)') 'Max bmlt (m/y), i, j     ',   &
            max_bmlt_global*scyr, imax_global, jmax_global
       call write_log(trim(message), type = GM_DIAGNOSTIC)
    endif

    ! max basal melt rate for grounded ice
    imax = 0
    jmax = 0
    max_bmlt_ground = unphys_val

    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%basal_melt%bmlt_ground(i,j) > max_bmlt_ground) then
             max_bmlt_ground = model%basal_melt%bmlt_ground(i,j)
             imax = i
             jmax = j
          endif
       enddo
    enddo

    call parallel_reduce_maxloc(xin=max_bmlt_ground, xout=max_bmlt_ground_global, xprocout=procnum)
    call parallel_globalindex(imax, jmax, imax_global, jmax_global, parallel)
    call broadcast(imax_global, proc = procnum)
    call broadcast(jmax_global, proc = procnum)

    ! Write to diagnostics only if nonzero

    if (abs(max_bmlt_global*scyr) > eps) then
       write(message,'(a25,f24.16,2i6)') 'Max bmlt grnd (m/y), i, j',   &
            max_bmlt_ground_global*scyr, imax_global, jmax_global
       call write_log(trim(message), type = GM_DIAGNOSTIC)
    endif

    ! max surface speed
    imax = 0
    jmax = 0
    max_spd_sfc = unphys_val

    do j = lhalo+1, velo_ns_ubound
       do i = lhalo+1, velo_ew_ubound
          spd = sqrt(model%velocity%uvel(1,i,j)**2   &
                   + model%velocity%vvel(1,i,j)**2)
          if (model%geomderv%stagthck(i,j) > minthck .and. spd > max_spd_sfc) then
             max_spd_sfc = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo

    call parallel_reduce_maxloc(xin=max_spd_sfc, xout=max_spd_sfc_global, xprocout=procnum)
    call parallel_globalindex(imax, jmax, imax_global, jmax_global, parallel)
    call broadcast(imax_global, proc = procnum)
    call broadcast(jmax_global, proc = procnum)

    write(message,'(a25,f24.16,2i6)') 'Max sfc spd (m/y), i, j  ',   &
                    max_spd_sfc_global*scyr, imax_global, jmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max basal speed

    imax = 0
    jmax = 0
    max_spd_bas = unphys_val
    do j = lhalo+1, velo_ns_ubound
       do i = lhalo+1, velo_ew_ubound
          spd = sqrt(model%velocity%uvel(upn,i,j)**2   &
                   + model%velocity%vvel(upn,i,j)**2)
          if (model%geomderv%stagthck(i,j) > minthck  .and. spd > max_spd_bas) then
             max_spd_bas = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo

    call parallel_reduce_maxloc(xin=max_spd_bas, xout=max_spd_bas_global, xprocout=procnum)
    call parallel_globalindex(imax, jmax, imax_global, jmax_global, parallel)
    call broadcast(imax_global, proc = procnum)
    call broadcast(jmax_global, proc = procnum)
    write(message,'(a25,f24.16,2i6)') 'Max base spd (m/y), i, j ',   &
                    max_spd_bas_global*scyr, imax_global, jmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    !-----------------------------------------------------------------
    ! Optionally, compute the rms_error of model thickness and surface speed
    ! compared to observations.
    ! TODO - Instead of ice_mask, use a mask that also includes cells
    !        which are ice-free in the model but ice-covered in observations
    !-----------------------------------------------------------------

    if (.not.parallel_is_zero(model%geometry%usrf_obs)) then

       call glissade_usrf_to_thck(&
            model%geometry%usrf_obs,      &
            model%geometry%topg,          &
            model%climate%eus,            &
            thck_obs)

       call glissade_rms_error(&
            ewn,            nsn,          &
            ice_mask,                     &
            parallel,                     &
            model%geometry%thck,          &  ! m
            thck_obs,                     &  ! m
            rmse_thck)

       write(message,'(a25,f24.16)') 'rms error, thickness (m) ', rmse_thck
       call write_log(trim(message), type = GM_DIAGNOSTIC)

    endif

    if (.not.parallel_is_zero(model%velocity%velo_sfc_obs)) then

       velo_sfc = sqrt(model%velocity%uvel(1,:,:)**2   &
                     + model%velocity%vvel(1,:,:)**2)

       call glissade_rms_error(&
            ewn-1,          nsn-1,              &
            stag_ice_mask,                      &
            parallel,                           &
            velo_sfc * scyr,                    &  ! m/yr
            model%velocity%velo_sfc_obs * scyr, &  ! m/yr
            rmse_velo)

       write(message,'(a25,f24.16)') 'rms error, sfc spd (m/y) ', rmse_velo
       call write_log(trim(message), type = GM_DIAGNOSTIC)

    endif

    ! local diagnostics

    ! initialize to unphysical negative values
    usrf_diag     = unphys_val
    thck_diag     = unphys_val
    topg_diag     = unphys_val
    relx_diag     = unphys_val
    load_diag     = unphys_val
    artm_diag     = unphys_val
    acab_diag     = unphys_val
    bmlt_diag     = unphys_val
    bwat_diag     = unphys_val
    bheatflx_diag = unphys_val
    top_age_diag  = unphys_val
    bot_age_diag  = unphys_val
    temp_diag(:)  = unphys_val
    spd_diag(:)   = unphys_val
    age_diag(:)   = unphys_val
    lithtemp_diag(:) = unphys_val    

    ! Set local diagnostic values, and communicate them to main_task
       
    if (model%numerics%idiag_local >= 1 .and. model%numerics%idiag_local <= ewn  &
                                        .and.                                    &
        model%numerics%jdiag_local >= 1 .and. model%numerics%jdiag_local <= nsn) then

       if (this_rank == model%numerics%rdiag_local) then

          i = model%numerics%idiag_local
          j = model%numerics%jdiag_local
          usrf_diag = model%geometry%usrf(i,j)
          thck_diag = model%geometry%thck(i,j)
          topg_diag = model%geometry%topg(i,j)
          if (model%options%isostasy == ISOSTASY_COMPUTE) then
             relx_diag = model%isostasy%relx(i,j)
             load_diag = model%isostasy%load(i,j)
          endif
          artm_diag = model%climate%artm_corrected(i,j)  ! artm_corrected = artm + artm_anomaly
          acab_diag = model%climate%acab_applied(i,j) * scyr
          bmlt_diag = model%basal_melt%bmlt_applied(i,j) * scyr
          if (model%options%which_ho_bwat == HO_BWAT_FLUX_ROUTING) then
             bwat_diag = model%basal_hydro%bwat_diag(i,j)
          else
             bwat_diag = model%basal_hydro%bwat(i,j)
          endif
          bheatflx_diag = model%temper%bheatflx(i,j)
          if (associated(model%geometry%ice_age)) then
             top_age_diag = model%geometry%ice_age(1,i,j)       ! age of top ice layer
             bot_age_diag = model%geometry%ice_age(upn-1,i,j)   ! age of bottom ice layer
          endif
          temp_diag(:) = model%temper%temp(1:upn,i,j)          
          spd_diag(:) = sqrt(model%velocity%uvel(1:upn,i,j)**2   &
                           + model%velocity%vvel(1:upn,i,j)**2) * scyr
          if (model%options%which_ho_ice_age == HO_ICE_AGE_COMPUTE) &
               age_diag(:) = model%geometry%ice_age(:,i,j)/scyr
          if (model%options%gthf == GTHF_COMPUTE) &
               lithtemp_diag(:) = model%lithot%temp(i,j,:)
       endif

       usrf_diag = parallel_reduce_max(usrf_diag)
       thck_diag = parallel_reduce_max(thck_diag)
       topg_diag = parallel_reduce_max(topg_diag)
       if (model%options%isostasy == ISOSTASY_COMPUTE) then
          relx_diag = parallel_reduce_max(relx_diag)
          load_diag = parallel_reduce_max(load_diag)
       endif
       artm_diag = parallel_reduce_max(artm_diag)
       acab_diag = parallel_reduce_max(acab_diag)
       bmlt_diag = parallel_reduce_max(bmlt_diag)
       bwat_diag = parallel_reduce_max(bwat_diag)
       bheatflx_diag = parallel_reduce_max(bheatflx_diag)
       if (associated(model%geometry%ice_age)) then
          top_age_diag = parallel_reduce_max(top_age_diag)
          bot_age_diag = parallel_reduce_max(bot_age_diag)
       endif
       do k = 1, upn
          temp_diag(k) = parallel_reduce_max(temp_diag(k))
          spd_diag(k)  = parallel_reduce_max(spd_diag(k))
       enddo
       do k = 1, upn-1
          age_diag(k)  = parallel_reduce_max(age_diag(k))
       enddo
       do k = 1, nlith
          lithtemp_diag(k) = parallel_reduce_max(lithtemp_diag(k))
       enddo

       call write_log(' ')
       write(message,'(a39,2i6)')  &
            'Grid point diagnostics: (i,j) =', model%numerics%idiag, &
                                               model%numerics%jdiag
       call write_log(trim(message), type = GM_DIAGNOSTIC)
       write(message,'(a39,3i6)')  &
            'Local (i,j,rank) =             ', model%numerics%idiag_local, &
                                               model%numerics%jdiag_local, &
                                               model%numerics%rdiag_local
       call write_log(trim(message), type = GM_DIAGNOSTIC)
       call write_log(' ')
 
       write(message,'(a25,f24.16)') 'Upper surface (m)        ', usrf_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Thickness (m)            ', thck_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Bedrock topo (m)         ', topg_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       if (model%options%isostasy == ISOSTASY_COMPUTE) then
          write(message,'(a25,f24.16)') 'Relaxed bedrock (m)   ', relx_diag
          call write_log(trim(message), type = GM_DIAGNOSTIC)
          write(message,'(a25,f24.16)') 'Load deflection (m)   ', load_diag
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       endif

       write(message,'(a25,f24.16)') 'Sfc mass balance (m/yr)  ', acab_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Basal mass balance (m/yr)', -bmlt_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Basal water depth (m)    ', bwat_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Basal heat flux (W/m^2)  ', bheatflx_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       ! Commented out, since we are writing the age in every column below
!       if (model%options%which_ho_ice_age == HO_ICE_AGE_COMPUTE) then
!          write(message,'(a25,f24.16)') 'Age of top layer (yr)    ', top_age_diag/scyr
!          call write_log(trim(message), type = GM_DIAGNOSTIC)
!          write(message,'(a25,f24.16)') 'Age of bottom layer (yr) ', bot_age_diag/scyr
!          call write_log(trim(message), type = GM_DIAGNOSTIC)
!       endif

       ! Vertical profile of ice speed and temperature

       call write_log(' ')
       write(message,'(a74)') ' Sigma       Ice speed (m/yr)      Ice temperature (C)        Ice age (yr)'
       call write_log(trim(message), type = GM_DIAGNOSTIC)
 
       if (size(model%temper%temp,1) == upn+1) then   ! temperatures staggered in vertical
                                                      ! (at layer midpoints)

           ! upper surface 
           write (message,'(f6.3,2f24.16)') model%numerics%sigma(1), spd_diag(1), artm_diag
           call write_log(trim(message), type = GM_DIAGNOSTIC)

           ! internal
           do k = 1, upn-1

              ! speed at top of layer
              if (k > 1) then
                 write (message,'(f6.3,f24.16)') model%numerics%sigma(k), spd_diag(k)
                 call write_log(trim(message), type = GM_DIAGNOSTIC)
              endif

              ! temp (and optionally age) at layer midpoint
              if (model%options%which_ho_ice_age == HO_ICE_AGE_COMPUTE) then
                 write (message,'(f6.3,24x,f24.16,f18.6)') model%numerics%stagsigma(k), temp_diag(k), age_diag(k)
                 call write_log(trim(message), type = GM_DIAGNOSTIC)
              else
                 write (message,'(f6.3,24x,f24.16)') model%numerics%stagsigma(k), temp_diag(k)
                 call write_log(trim(message), type = GM_DIAGNOSTIC)
              endif

           enddo

           ! lower surface
           write (message,'(f6.3,2f24.16)') model%numerics%sigma(upn), spd_diag(upn), temp_diag(upn)
           call write_log(trim(message), type = GM_DIAGNOSTIC)
           
       else    ! temperatures unstaggered in vertical (at layer interfaces)
 
           do k = 1, upn
             write (message,'(f6.3,2f24.16)') model%numerics%sigma(k), spd_diag(k), temp_diag(k)
             call write_log(trim(message), type = GM_DIAGNOSTIC)
          enddo

       endif  ! temps staggered

       ! Vertical profile of upper lithosphere temperature

       if (model%options%gthf == GTHF_COMPUTE) then

          call write_log(' ')
          write(message,'(a41)') '  Level (m)          Lithosphere temp (C)'
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          level = 0.d0
          do k = 1, nlith
             level = level + model%lithot%deltaz(nlith)
             write (message,'(f10.0,6x,f24.16)') level, lithtemp_diag(k)
             call write_log(trim(message), type = GM_DIAGNOSTIC)
          enddo

       endif  ! gthf_compute

    endif     ! idiag_local and jdiag_local in bounds

    call write_log(' ')

    ! glacier diagnostics

    if (model%options%enable_glaciers) then

       tot_glc_area = 0.0d0
       tot_glc_volume = 0.0d0
       tot_glc_area_init = 0.0d0
       tot_glc_volume_init = 0.0d0
       tot_glc_area_init_extent = 0.0d0
       tot_glc_volume_init_extent = 0.0d0
       count_area = 0
       count_volume = 0

       do ng = 1, model%glacier%nglacier
          tot_glc_area = tot_glc_area + model%glacier%area(ng)
          tot_glc_volume = tot_glc_volume + model%glacier%volume(ng)
          tot_glc_area_init = tot_glc_area_init + model%glacier%area_init(ng)
          tot_glc_volume_init = tot_glc_volume_init + model%glacier%volume_init(ng)
          tot_glc_area_init_extent = tot_glc_area_init_extent + model%glacier%area_init_extent(ng)
          tot_glc_volume_init_extent = tot_glc_volume_init_extent + model%glacier%volume_init_extent(ng)
          if (model%glacier%area(ng) > eps) then
             count_area = count_area + 1
          endif
          if (model%glacier%volume(ng) > eps) then
             count_volume = count_volume + 1
          endif
       enddo

       ! Copy selected scalars into the derived type
       model%glacier%total_area = tot_glc_area
       model%glacier%total_volume = tot_glc_volume
       model%glacier%nglacier_active = count_area

       ! Write some total glacier diagnostics

       write(message,'(a25)') 'Glacier diagnostics: '
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       call write_log(' ')

       write(message,'(a35,i14)')   'Number of glaciers                 ', &
            model%glacier%nglacier
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a35,i14)')   'Glaciers with nonzero area         ', &
            count_area
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a35,i14)')   'Glaciers with nonzero volume       ', &
            count_volume
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a35,f14.6)') 'Total glacier area_init (km^2)     ', &
            tot_glc_area_init / 1.0d6
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a35,f14.6)') 'Total glacier area (km^2)          ', &
            tot_glc_area / 1.0d6
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a35,f14.6)') 'Total area_init_extent (km^2)      ', &
            tot_glc_area_init_extent / 1.0d6
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a35,f14.6)') 'Total glacier volume_init (km^3)   ', &
            tot_glc_volume_init / 1.0d9
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a35,f14.6)') 'Total glacier volume (km^3)        ', &
            tot_glc_volume / 1.0d9
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a35,f14.6)') 'Total volume_init_extent (km^3)    ', &
            tot_glc_volume_init_extent / 1.0d9
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       if (model%glacier%set_powerlaw_c == GLACIER_POWERLAW_C_INVERSION) then

          ! diagnostics related to thickness inversion

          tot_glc_area_target = 0.0d0
          tot_glc_volume_target = 0.0d0
          do ng = 1, model%glacier%nglacier
             tot_glc_area_target = tot_glc_area_target + model%glacier%area_target(ng)
             tot_glc_volume_target = tot_glc_volume_target + model%glacier%volume_target(ng)
          enddo

          ! Compute the root-mean-square error of (thck - thck_target), including cells
          !  with cism_glacier_id > 0 or cism_glacier_id_init > 0
          where (model%glacier%cism_glacier_id_init > 0 .or. model%glacier%cism_glacier_id > 0)
             glacier_ice_mask = 1
          elsewhere
             glacier_ice_mask = 0
          endwhere

          call glissade_rms_error(&
               ewn,            nsn,          &
               glacier_ice_mask,             &
               parallel,                     &
               model%geometry%thck,          &
               model%glacier%thck_target,    &
               glc_rmse_thck)

          ! Repeat, including only cells within the initial glacier extent

          where (model%glacier%cism_glacier_id_init > 0)
             glacier_ice_mask = 1.0d0
          elsewhere
             glacier_ice_mask = 0.0d0
          endwhere

          call glissade_rms_error(&
               ewn,            nsn,          &
               glacier_ice_mask,             &
               parallel,                     &
               model%geometry%thck,          &
               model%glacier%thck_target,    &
               glc_rmse_thck_init_extent)

          write(message,'(a35,f14.6)') 'Total area target (km^2)           ', &
               tot_glc_area_target / 1.0d6
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'Total volume target (km^2)         ', &
               tot_glc_volume_target / 1.0d9
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'rms error, thck - thck_target (m)  ', &
               glc_rmse_thck
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'rms error over init extent (m)     ', &
               glc_rmse_thck_init_extent
          call write_log(trim(message), type = GM_DIAGNOSTIC)

       endif  ! set_powerlaw_c

       call write_log(' ')

       ! Write output related to the diagnostic glacier

       ng = model%glacier%ngdiag

       if (ng > 0) then

          write(message,'(a35,i14)') 'Diagnostic glacier index (RGI)     ', &
               model%glacier%cism_to_rgi_glacier_id(ng)
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,i14)') 'Diagnostic glacier index (CISM)    ', ng
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'Glacier area_init (km^2)           ', &
               model%glacier%area_init(ng) / 1.0d6
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'Glacier area (km^2)                ', &
               model%glacier%area(ng) / 1.0d6
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'Glacier area_init_extent (km^2)    ', &
               model%glacier%area_init_extent(ng) / 1.0d6
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'Glacier volume (km^3)              ', &
               model%glacier%volume(ng) / 1.0d9
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'Glacier volume_init (km^3)         ', &
               model%glacier%volume_init(ng) / 1.0d9
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'Glacier volume_init_extent (km^3)  ', &
               model%glacier%volume_init_extent(ng) / 1.0d9
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'mu_star (mm/yr w.e./deg C)         ', &
               model%glacier%mu_star(ng)
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'alpha_snow                         ', &
               model%glacier%alpha_snow(ng)
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a35,f14.6)') 'beta_artm (deg C)                  ', &
               model%glacier%beta_artm(ng)
          call write_log(trim(message), type = GM_DIAGNOSTIC)

       else   ! glacier ID = 0

          write(message,'(a35,i14)') 'Diagnostic glacier index (CISM)    ', ng
          call write_log(trim(message), type = GM_DIAGNOSTIC)

       endif

       call write_log(' ')

    endif  ! enable_glaciers

  end subroutine glissade_write_diag

!****************************************************************************

  subroutine glissade_mass_balance_diagnostics(model)

    ! Compute mass balance diagnostics associated with five processes:
    ! (1) surface mass balance, (2) basal mass balance, (3) calving,
    ! (4) lateral melt, and (5) other kinds of ice removal.
    ! Ice cap removal, if applied, falls under (5).

    use glimmer_physcon, only: rhoi, rhow

    type(glissade_global_type), intent(inout) :: model   ! model instance

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
    model%mass_flux%latmelt_flux(:,:) = rhoi * (-model%lateral_melt%latmelt_thck(:,:)) / model%numerics%dt
    model%mass_flux%removal_flux(:,:) = rhoi * (-model%geometry%removal_thck(:,:)) / model%numerics%dt

    ! rates of calving, lateral melt and ice removal (m/yr ice; positive for ice loss)
    model%calving%calving_rate(:,:) = model%calving%calving_thck(:,:) / (model%numerics%dt/scyr)
    model%lateral_melt%latmelt_rate(:,:) = model%lateral_melt%latmelt_thck(:,:) / (model%numerics%dt/scyr)
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

    type(glissade_tensor_type), intent(in) :: &
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

    type(glissade_tensor_type), intent(inout) :: &
         strain_rate               ! 3D strain rate tensor
                                   ! intent(out) if computed from tau and efvs

    real(dp), dimension(nx,ny), intent(out) :: &
         eps_eigen1, eps_eigen2    ! eigenvalues of 2D horizontal stress tensor (1/s)

    type(glissade_tensor_type), intent(in), optional :: &
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

    type(glissade_global_type), intent(inout) :: model   ! model instance

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
