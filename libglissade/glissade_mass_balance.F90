!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_mass_balance.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains drivers for handling the surface and basal mass balance.
!
! Author: William Lipscomb
!         NSF National Center for Atmospheric Research
!         Climate and Global Dynamics Laboratory
!         Boulder, CO 80303
!         USA
!         <lipscomb@ucar.edu>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
  module glissade_mass_balance

    use glimmer_global, only: dp
    use glimmer_log
    use glimmer_utils, only: point_diag
    use glide_types
    use cism_parallel, only: this_rank, main_task, nhalo, lhalo, uhalo, &
         parallel_halo, parallel_reduce_max, parallel_reduce_sum, parallel_globalindex

    implicit none
    save
    private

    public :: glissade_mass_balance_init, glissade_prepare_climate_forcing,  &
         glissade_apply_smb, glissade_add_2d_anomaly, glissade_add_3d_anomaly
    public :: verbose_smb

!!    logical, parameter :: verbose_smb = .false.
    logical, parameter :: verbose_smb = .true.

    logical, parameter ::     &
         conservation_check = .true. ! if true, check global conservation

!=======================================================================

  contains

!=======================================================================

  subroutine glissade_mass_balance_init(model)

    ! Initialize some fields related to the surface mass balance

    use glimmer_paramets, only: eps11
    use glimmer_physcon, only: rhow, rhoi, scyr

    ! input/output arguments

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! local variables

    real(dp) :: local_maxval, global_maxval
    character(len=100) :: message
      
    ! Initialize acab, if SMB (with different units) was read in
    if (model%options%smb_input == SMB_INPUT_MMYR_WE) then
       ! Convert units from mm/yr w.e. to m/s ice
       model%climate%acab(:,:) = (model%climate%smb(:,:) * (rhow/rhoi)/1000.d0) / scyr
       !WHL - debug
       if (main_task) write(6,*) 'Setting acab, m/yr ice'
    endif

    ! Initialize artm_corrected.  This is equal to artm, plus any prescribed temperature anomaly.
    model%climate%artm_corrected(:,:) = model%climate%artm(:,:)
    
    if (model%options%enable_artm_anomaly) then
       ! Check whether artm_anomaly was read from an external file.
       ! If so, then use this field as the anomaly.
       ! If not, then set artm_anomaly = artm_anomaly_constant everywhere.
       ! Note: The artm_anomaly field does not change during the run,
       !       but it is possible to ramp up the anomaly using artm_anomaly_timescale.

       local_maxval = maxval(abs(model%climate%artm_anomaly))
       global_maxval = parallel_reduce_max(local_maxval)
       if (global_maxval < eps11) then
          model%climate%artm_anomaly = model%climate%artm_anomaly_const
          write(message,*) &
               'Setting artm_anomaly = constant value (degC):', model%climate%artm_anomaly_const
          call write_log(trim(message))
       else
          if (model%options%is_restart == NO_RESTART) then
             call write_log('Setting artm_anomaly from external file')
          endif
       endif
    endif
    !TODO - Write a short utility function to compute global_maxval of any field.
    !TODO - Repeat for snow and precip anomalies

    ! If acab is to be overwritten for some cells, then set overwrite_acab_mask = 1 for these cells.
    ! We can overwrite the input acab with a fixed value (typically negative) where
    ! (1) the input acab = 0 at initialization, or
    ! (2) the input thck <= overwrite_acab_minthck at initialization
    ! Note: This option is designed for standalone runs, and should be used with caution for coupled runs.
    !       On restart, overwrite_acab_mask is read from the restart file.
    
    if (model%climate%overwrite_acab_value /= 0 .and. model%options%is_restart == NO_RESTART) then

       call set_overwrite_acab_mask(&
            model%options%overwrite_acab,          &
            model%climate%acab,                    &
            model%geometry%thck,                   &
            model%climate%overwrite_acab_minthck,  &
            model%climate%overwrite_acab_mask)

    endif
  
  end subroutine glissade_mass_balance_init

!=======================================================================

  subroutine glissade_prepare_climate_forcing(model)

    ! Given the climate input fields (some combination of artm, acab, snow and precip)
    !  and the SMB input options (2D, 2D with vertical gradient, or 3D), compute the current SMB,
    !  downscaling to the current surface elevation if needed.
    ! Also add anomaly terms, if needed, and do any required unit conversions.

    use glimmer_physcon, only: rhow, rhoi, scyr
    use glissade_grid_operators, only: glissade_vertical_interpolate
    use glissade_masks, only: glissade_extend_mask
    use cism_parallel, only: parallel_is_nonzero

    ! input/output arguments

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! local variables

    integer, dimension(model%general%ewn, model%general%nsn) ::   &
       extended_ice_sheet_mask   ! extension of ice_sheet_mask to include neighbor cells

    integer :: itest, jtest, rtest      ! coordinates of diagnostic cell
    integer :: i, j, k
    integer :: ewn, nsn, upn
    integer :: nlev_smb                 ! number of levels at which the 3D SMB is computed

    type(parallel_type) :: parallel     ! info for parallel communication

    character(len=100) :: message

    ! initialize

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn
    nlev_smb = model%climate%nlev_smb

    parallel = model%parallel

    !-------------------------------------------------------
    !WHL - Notes for Ishraque
    ! This could be a good place to call your paleo interpolation subroutine.
    ! For now, you could trigger the call based on model%climate%input_function == SMB_INPUT_FUNCTION_PDD,
    !  although in the long run the paleo interpolation should be a separate option.
    !
    ! Inputs:
    !   climate index, a scalar which will have been read from a forcing file and updated based on the current time
    !   artm_ref for each of the three climates, based on the current time of year
    !
    ! Outputs:
    !   model%climate%artm_ref for the current time, obtained by paleo interpolation
    !   model%climate%artm_gradz if artm_input_function = ARTM_INPUT_FUNCTION_XY_GRADZ = 1, obtained by paleo interpolation
    !     (not needed if artm_input_function = ARTM_INPUT_FUNCTION_XY_LAPSE = 3)
    !   model%climate%precip for the current time, based on artm_ref and Clausius-Clapeyron
    !
    ! Maybe there are other fields I've forgotten?
    ! I think you will read model%climate%usrf_ref at initialization.
    !
    ! The parameters degree_factor, tmlt, snow_threshold_min, snow_threshold_max and t_lapse will be set in the config file.
    !
    ! Downscaling of artm_ref to artm (at the ice surface) happens below, followed by the SMB calculation.
    !-------------------------------------------------------

    ! Downscale artm to the current surface elevation if needed.
    ! The downscaling options are:
    ! (0) artm(x,y); no dependence on surface elevation
    ! (1) artm(x,y) + d(artm)/dz(x,y) * dz; artm depends on input field at reference elevation, plus vertical correction
    ! (2) artm(x,y,z); artm obtained by linear interpolation between values prescribed at adjacent vertical levels
    ! (3) artm(x,y) adjusted with a uniform lapse rate
    ! For options (1) - (3), the elevation-dependent artm is computed here.

    if (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XY_GRADZ) then

       ! compute artm by a lapse-rate correction to the reference value
       model%climate%artm(:,:) = model%climate%artm_ref(:,:) + &
            (model%geometry%usrf(:,:) - model%climate%usrf_ref(:,:)) * model%climate%artm_gradz(:,:)

    elseif (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XYZ) then

       ! Note: With linear_extrapolate_in = T, the values outside the range are obtained by linear extrapolation
       !        from the top two or bottom two values.
       !       For temperature, which varies roughly linearly with elevation, this is more accurate
       !        than simply extending the top and bottom values.
       !       This call includes a halo update.

       call glissade_vertical_interpolate(ewn,                   nsn,                       &
                                          nlev_smb,              model%climate%smb_levels,  &
                                          model%geometry%usrf,                              &
                                          model%climate%artm_3d,                            &
                                          model%climate%artm,                               &
                                          linear_extrapolate_in = .true.)

    elseif (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XY_LAPSE) then

       ! compute artm by a lapse-rate correction to artm_ref
       ! T_lapse is defined as positive for T decreasing with height

       model%climate%artm(:,:) = model%climate%artm_ref(:,:) - &
            (model%geometry%usrf(:,:) - model%climate%usrf_ref(:,:)) * model%climate%t_lapse

    endif   ! artm_input_function

    call parallel_halo(model%climate%artm, parallel)

    ! Optionally, add an anomaly to the surface air temperature
    ! Typically, artm_corrected = artm, but sometimes (e.g., for ISMIP6 forcing experiments),
    !  it includes a time-dependent anomaly.
    ! Note that artm itself does not change in time, unless it is elevation-dependent.

    model%climate%artm_corrected(:,:) = model%climate%artm(:,:)

    if (model%options%enable_artm_anomaly) then

       call glissade_add_2d_anomaly(&
            model%climate%artm_corrected,          &   ! degC
            model%climate%artm_anomaly,            &   ! degC
            model%climate%artm_anomaly_tstart,     &   ! yr
            model%climate%artm_anomaly_timescale,  &   ! yr
            model%numerics%time)                       ! yr

    endif

    ! Similar calculations for snow and precip anomalies
    ! Note: These variables are currently used only to compute glacier SMB.
    !       There are assumed to have the same timescale as artm_anomaly.
    ! TODO: Define a single anomaly timescale for all anomaly forcing?

    model%climate%snow_corrected(:,:) = model%climate%snow(:,:)

    if (model%options%enable_snow_anomaly) then

       call glissade_add_2d_anomaly(&
            model%climate%snow_corrected,          &   ! mm/yr w.e.
            model%climate%snow_anomaly,            &   ! mm/yr w.e.
            model%climate%artm_anomaly_tstart,     &   ! yr
            model%climate%artm_anomaly_timescale,  &   ! yr
            model%numerics%time)                       ! yr

    endif

    model%climate%precip_corrected(:,:) = model%climate%precip(:,:)

    if (model%options%enable_precip_anomaly) then

       call glissade_add_2d_anomaly(&
            model%climate%precip_corrected,        &   ! mm/yr w.e.
            model%climate%precip_anomaly,          &   ! mm/yr w.e.
            model%climate%artm_anomaly_tstart,     &   ! yr
            model%climate%artm_anomaly_timescale,  &   ! yr
            model%numerics%time)                       ! yr

    endif

    if (verbose_smb .and. this_rank==rtest) then
       if (model%options%enable_artm_anomaly) then
          i = itest
          j = jtest
          write(6,*) 'rank, i, j, time, anomaly timescale (yr):', &
               this_rank, i, j, model%numerics%time, model%climate%artm_anomaly_timescale
          write(6,*) '   artm, artm anomaly, corrected artm (deg C):', model%climate%artm(i,j), &
               model%climate%artm_anomaly(i,j), model%climate%artm_corrected(i,j)
          if (model%options%enable_snow_anomaly) then
             write(6,*) '   snow, snow anomaly, corrected snow (mm/yr):', model%climate%snow(i,j), &
                  model%climate%snow_anomaly(i,j), model%climate%snow_corrected(i,j)
          endif
          if (model%options%enable_precip_anomaly) then
             write(6,*) '   prcp, prcp anomaly, corrected prcp (mm/yr):', model%climate%precip(i,j), &
                  model%climate%precip_anomaly(i,j), model%climate%precip_corrected(i,j)
          endif
       endif   ! enable_artm_anomaly
    endif   ! verbose

    !-------------------------------------------------------------------------
    ! Some notes on SMB fields and units:
    !
    ! Most of the SMB calculations in Glissade are based on model%climate%acab,
    !  defined as the net forcing (accumulation minus ablation) at the upper ice surface.
    ! Not all of this SMB is necessarily applied; for instance, a negative SMB cannot
    !  be applied to ice-free cells. The applied SMB is stored in model%climate%acab_applied.
    !
    ! Both acab and acab_applied have units of m/yr ice in I/O, and m/s ice in the code.
    ! They must be multiplied by rhoi/rhow to convert to m/s w.e.
    !
    ! If the SMB is passed from CESM, then it arrives in Glad with units of kg/m^2/s w.e.,
    !  and Glad does the conversions needed to compute model%climate%acab in units of m/s ice.
    !
    ! For standalone ice runs, the input SMB is often supplied in units of kg/m^2/yr w.e.,
    !  or equivalently mm w.e./yr. Sometimmes it is supplied at a reference elevation
    !  and must be downscaled to the ice surface, possibly with a vertical gradient.
    !  Or it can be supplied as a 3D field at multiple levels, and then vertically interpolated.
    !  In other runs, the SMB may be suppied as a baseline field to which an anomaly is added.
    !
    ! Subroutine glissade_apply_smb handles these various cases, if applicable.
    ! It assumes that smb and related fields (smb_anomaly, etc.) have units of mm w.e./yr.
    ! After doing the desired adjustments, we have a field called smb_corrected.
    ! This field is multiplied by (rhow/rhoi)/scyr to compute model%climate%acab,
    !  which is passed to the main mass balance driver.
    !-------------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Depending on the SMB input options, compute model%climate%acab at the ice surface.
    ! The options are:
    ! (0) SMB(x,y); no dependence on surface elevation
    ! (1) SMB(x,y) + dSMB/dz(x,y) * dz; SMB depends on input field at reference elevation, plus vertical correction
    ! (2) SMB(x,y,z); SMB obtained by linear interpolation between values prescribed at adjacent vertical levels
    ! (3) SMB obtained from precip and artm using a positive-degree scheme
    !
    ! Options (1) and (2) require input fields with SMB units of mm/yr w.e. (SMB_INPUT_MMYR_WE)
    ! For these options, the elevation-dependent SMB is computed here.
    ! ------------------------------------------------------------------------

    if (model%options%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ) then

       ! downscale SMB to the local surface elevation
       model%climate%smb(:,:) = model%climate%smb_ref(:,:) + &
            (model%geometry%usrf(:,:) - model%climate%usrf_ref(:,:)) * model%climate%smb_gradz(:,:)

    elseif (model%options%smb_input_function == SMB_INPUT_FUNCTION_XYZ) then

       ! downscale SMB to the local surface elevation
       ! Note: With linear_extrapolate_in = F, the values at top and bottom levels are simply extended upward and downward.
       !       For SMB, this is safer than linear extrapolation (especially when extrapolating upward).

       call glissade_vertical_interpolate(&
            ewn,       nsn,                       &
            nlev_smb,  model%climate%smb_levels,  &
            model%geometry%usrf,                  &
            model%climate%smb_3d,                 &
            model%climate%smb,                    &
            linear_extrapolate_in = .false.)

    elseif  (model%options%smb_input_function == SMB_INPUT_FUNCTION_PDD) then

       ! Compute SMB using a simple PDD scheme:
       ! (1) Partition precip as rain or snow based on the downscaled artm
       ! (2) Compute ablation based on artm and a degree factor
       ! Assume that artm has already been downscaled, if needed, based on artm_input_function.

       ! Note: This is similar to the SMB calculation for glaciers, but that calculation is done in the glacier module.
       ! TODO: Put the glacier values of snow_threshold_min and snow_threshold_max in the climate derived type.

       ! compute snow accumulation (mm/yr w.e.)
       where (model%climate%artm > model%climate%snow_threshold_max)
          model%climate%snow = 0.0d0   ! all precip falls as rain
       elsewhere (model%climate%artm < model%climate%snow_threshold_min)
          model%climate%snow = model%climate%precip   ! all precip falls as snow
       elsewhere (model%climate%artm > model%climate%snow_threshold_min)
          model%climate%snow = model%climate%precip * (model%climate%snow_threshold_max - model%climate%artm)  &
               / (model%climate%snow_threshold_max - model%climate%snow_threshold_min)
       endwhere

       ! compute ablation (mm/yr w.e.)
       ! Note: degree_factor has units of mm/yr w.e./degC to be consistent with other mass-balance variables.
       !       It is like mu_star for glaciers.
       model%climate%ablation = model%climate%degree_factor * max(model%climate%artm - model%climate%tmlt, 0.0d0)

       ! compute smb (mm/yr w.e.)
       model%climate%smb = model%climate%snow - model%climate%ablation

       ! set smb = 0 for open ocean
       where (model%geometry%thck == 0.0d0 .and. (model%geometry%topg - model%climate%eus) < 0.0d0)
          model%climate%smb = 0.0d0
       endwhere

    endif   ! smb_input_function

    ! For the non-default smb_input_function options, make sure that model%climate%smb is nonzero somewhere; else abort.
    ! For the default option, do not abort, since idealized tests often have a zero SMB.

    call parallel_halo(model%climate%smb, parallel)

    if (model%options%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ .or. &
        model%options%smb_input_function == SMB_INPUT_FUNCTION_XYZ .or. &
        model%options%smb_input_function == SMB_INPUT_FUNCTION_PDD) then
       if (parallel_is_nonzero(model%climate%smb)) then
          ! all is well
       else
          write(message,*) 'Error: smb = 0 everywhere with smb_input_function =', model%options%smb_input_function
          call write_log(trim(message), GM_FATAL)
       endif
    endif

    ! optional diagnostics
    if (verbose_smb) then

       if (this_rank == rtest) then
          write(6,*) 'Computing runtime smb with smb_input_function =', model%options%smb_input_function
       endif
       call point_diag(model%geometry%usrf, 'usrf (m)', itest, jtest, rtest, 7, 7)

       if (model%options%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ) then
          call point_diag(model%geometry%usrf - model%climate%usrf_ref, 'usrf - usrf_ref (m)', &
               itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%smb_ref, 'reference smb (mm/yr)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%smb_gradz, 'smb_gradz (mm/yr per m)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%smb, 'downscaled smb (mm/yr)', itest, jtest, rtest, 7, 7)
       elseif (model%options%smb_input_function == SMB_INPUT_FUNCTION_XYZ) then
          k = model%climate%nlev_smb/2 + 1  ! arbitrary k
          if (this_rank == rtest) write(6,*) 'Diagnostic level k, level (m) =', k, model%climate%smb_levels(k)
          call point_diag(model%climate%smb_3d(k,:,:), '3d smb (mm/yr)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%smb, 'downscaled smb (mm/yr)', itest, jtest, rtest, 7, 7)
       elseif (model%options%smb_input_function == SMB_INPUT_FUNCTION_PDD) then
          call point_diag(model%climate%artm, 'artm (deg C)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%precip, 'precip (mm/yr)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%snow, 'snow (mm/yr)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%ablation, 'ablation (mm/yr)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%smb,'smb (mm/yr)', itest, jtest, rtest, 7, 7)
       endif  ! smb_input_function

       if (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XY_GRADZ) then
          call point_diag(model%climate%artm_ref, 'reference artm (deg C)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%artm_gradz*1000.d0, 'artm_gradz (deg C per km)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%artm, 'downscaled artm (deg C)', itest, jtest, rtest, 7, 7)
       elseif (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XYZ) then
          k = model%climate%nlev_smb/2 + 1  ! arbitrary k
          if (this_rank == rtest) write(6,*) 'Diagnostic level k, level (m) =', k, model%climate%smb_levels(k)
          call point_diag(model%climate%artm_3d(k,:,:), '3d artm (deg C)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%artm, 'downscaled artm (deg C)', itest, jtest, rtest, 7, 7)
       elseif (model%options%artm_input_function == ARTM_INPUT_FUNCTION_XY_LAPSE) then
          call point_diag(model%climate%artm_ref, 'reference artm (deg C)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%artm, 'downscaled artm (deg C)', itest, jtest, rtest, 7, 7)
       endif   ! artm_input_function

    endif  ! verbose_smb

    ! Compute a corrected smb field that includes any anomalies or correction factors.

    ! initialize
    model%climate%smb_corrected(:,:) = model%climate%smb(:,:)

    ! Optionally, multiply smb by a scalar adjustment factor
    if (model%climate%smb_factor /= 1.0d0) then
       model%climate%smb_corrected(:,:) = model%climate%smb_corrected(:,:) * model%climate%smb_factor
    endif

    ! Optionally, add an anomaly field
    ! Note: If enable_artm_anomaly = T, then artm_anomaly is computed above, before calling glissade_therm_driver.

    if (model%options%enable_smb_anomaly) then

       call parallel_halo(model%climate%smb_anomaly, parallel)

       call glissade_add_2d_anomaly(&
            model%climate%smb_corrected,           &   ! scaled model units
            model%climate%smb_anomaly,             &   ! mm/yr w.e.
            model%climate%smb_anomaly_tstart,      &   ! yr
            model%climate%smb_anomaly_timescale,   &   ! yr
            model%numerics%time)                       ! yr

       if (verbose_smb .and. this_rank==rtest) then
          i = itest
          j = jtest
          write(6,*) 'i, j, time, input smb, smb anomaly, corrected smb (mm/yr):', &
               i, j, model%numerics%time, model%climate%smb(i,j), model%climate%smb_anomaly(i,j), model%climate%smb_corrected(i,j)
       endif

    endif

    ! Note on glacier SMB calculations
    ! Glaciers require smb_input = SMB_INPUT_MMYR_WE and compute SMB as follows:
    ! * Temperature and snowfall are accumulated during each call to glissade_glacier_update.
    ! * The annual mean SMB is computed at the end of the year in units of mm/yr w.e.
    !   and copied to model%climate%smb,
    ! * The SMB is then applied uniformly during the following year.
    ! Thus, the only thing to do here is to convert SMB (mm/yr w.e.) to acab (m/s ice).

    ! Convert units from mm/yr w.e. to acab units of m/s ice, if needed.
    if (model%options%smb_input == SMB_INPUT_MMYR_WE) then
       model%climate%acab(:,:) = (model%climate%smb_corrected(:,:) * (rhow/rhoi)/1000.d0) / scyr
    endif

    ! Optionally, overwrite acab where overwrite_acab_mask = 1.
    ! Note: overwrite_acab_value has scaled model units for now
    if (model%options%overwrite_acab /= 0) then
       call overwrite_acab(&
            model%climate%overwrite_acab_mask,  &
            model%climate%overwrite_acab_value, &
            model%climate%acab)
    endif

    !TODO - Compute ice_sheet_mask here?
    ! Optionally, block ice sheet inception by allowing SMB to be nonzero
    !  only in the main ice sheet and in cells adjacent to the ice sheet.
    ! Note: The ice sheet mask is computed in the diagnostic solve at the end of the previous time step
    !       (and at initialization).  Since it includes all cells that were part of the ice sheet
    !       before transport, an extended version of the mask will include all cells that potentially
    !       are part of the ice sheet after transport.

    if (model%options%block_inception) then

       ! Extend the ice sheet mask to include nearest neighbors (both edge and corner neighbors).

       call glissade_extend_mask(&
            model%general%ewn,   model%general%nsn,     &
            model%geometry%ice_sheet_mask,              &
            extended_mask = extended_ice_sheet_mask)

       call parallel_halo(model%geometry%ice_sheet_mask, parallel)

       ! SMB is allowed to be positive only where ice_mask_smb = 1.
       ! Note: This code allows cells outside the ice sheet to melt.

       where (extended_ice_sheet_mask == 0)
          model%climate%acab = min(model%climate%acab, 0.0d0)
       endwhere

    endif

    ! Optionally, correct acab by adding (-dthck_dt_obs_basin) where ice is floating.
    ! (Note: model%climate%acab does not change.)
    ! During inversions for deltaT_ocn, this will generally force a positive ocean melt rate
    !  where the ice is thinning, preventing large negative values of deltaT_ocn during spin-up.
    ! When the correction is removed, the ice should melt and thin in agreement with observations.
    ! Algorithm:
    ! (1) For each basin, compute the average of dthck_dt_obs over floating ice.
    !     Include all floating cells in the average.
    ! (2) For all cells in each basin, set dthck_dt_obs_basin to this average.
    !     Limit so that dthck_dt_obs_basin <= 0.
    ! (3) At runtime, add (-dthck_dt_obs_basin) to acab for each floating cell.

    if (model%options%enable_acab_dthck_dt_correction) then

       if (verbose_smb) then
          call point_diag(model%climate%acab*scyr, 'acab (m/yr)', itest, jtest, rtest, 7, 7)
       endif

       where (model%geometry%f_ground_cell < 1.0d0 .and. model%geometry%dthck_dt_obs_basin < 0.0d0)
          ! floating ice is thinning in obs; apply a positive correction to acab
          ! Note: dthck_dt_obs_basin has units of m/yr; convert to m/s
          model%climate%acab = model%climate%acab &
               - (1.0d0 - model%geometry%f_ground_cell) * (model%geometry%dthck_dt_obs_basin/scyr)
       endwhere

       if (verbose_smb) then
          call point_diag(-model%geometry%dthck_dt_obs_basin, 'dthck_dt_obs correction (m/yr)', itest, jtest, rtest, 7, 7)
          call point_diag(model%climate%acab*scyr, 'new acab (m/yr)', itest, jtest, rtest, 7, 7)
       endif

    endif   ! enable_acab_dthck_dt_correction

  end subroutine glissade_prepare_climate_forcing

!=======================================================================

  subroutine glissade_apply_smb(model)

    ! Apply the SMB at the upper and lower surfaces, and recompute tracer values.

    use glimmer_paramets, only: eps11
    use glimmer_physcon, only: rhow, rhoi, scyr
    use glissade_masks, only: glissade_get_masks
    use glissade_calving, only: verbose_calving

    ! input/output arguments

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! local variables

    real(dp), dimension(model%general%ewn,model%general%nsn) ::   &
         bmlt                 ! = model%basal_melt%bmlt (m/s) if basal mass balance is included in continuity equation, else = 0

    ! masks
    integer, dimension(model%general%ewn, model%general%nsn) ::   &
         ice_mask,               & ! = 1 if thck > 0, else = 0
         floating_mask,          & ! = 1 where ice is present and floating, else = 0
         ocean_mask,             & ! = 1 if topg is below sea level and thck = 0, else = 0
         land_mask,              & ! = 1 if topg is at or above sea level, else = 0
         calving_front_mask        ! = 1 where ice is floating and borders an ocean cell, else = 0

    character(len=100) :: message

    integer :: itest, jtest, rtest      ! coordinates of diagnostic cell
    integer :: i, j, k
    integer :: ewn, nsn, upn
    integer :: nlev_smb                 ! number of levels at which the 3D SMB is computed

    type(parallel_type) :: parallel     ! info for parallel communication

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn
    nlev_smb = model%climate%nlev_smb

    parallel = model%parallel

    ! Convert bmlt to SI units (m/s)
    ! Note: bmlt is the sum of bmlt_ground (computed in glissade_thermal_solve) and bmlt_float
    !       (computed in glissade_bmlt_float_solve).
    ! Note: bmlt can be turned off by setting options%basal_mbal = BASAL_MBAL_NO_CONTINUITY

    if (model%options%basal_mbal == BASAL_MBAL_CONTINUITY) then    ! include bmlt in continuity equation
       bmlt(:,:) = model%basal_melt%bmlt(:,:)
    else                                                           ! do not include bmlt in continuity equation
       bmlt = 0.0d0
    endif

    ! ------------------------------------------------------------------------
    ! Get masks used for the mass balance calculation.
    ! Pass thklim = 0 to identify cells with thck > 0 (not thck > thklim).
    ! Use ocean_mask to identify ocean cells where positive acab should not be applied.
    ! Use thck_effective to compute a fractional area for calving_front cells.
    ! TODO - Is it correct to use the old value of f_ground_cell from the start of the time step?
    !        Note that this value is used to identify CF cells where the mass balance is corrected.
    ! TODO - Would it be better not to recompute the masks here, but to use the pre-transport masks?
    ! ------------------------------------------------------------------------

    !Note: If not using the subgrid CF, then we should recompute effective_areafrac before
    !       the call to mass_balance_driver. This allows us to remove small thin ice
    !       from marine cells with negative acab or positive bmlt.
    !      If using the subgrid CF, use the pre-transport masks. Thin ice in ocean cells
    !       will be removed with different logic based on protected_mask.

    if (model%options%which_ho_calving_front == HO_CALVING_FRONT_NO_SUBGRID) then

       call glissade_get_masks(&
            ewn,              nsn,              &
            parallel,                           &
            model%geometry%thck,                &   ! m
            model%geometry%topg,                &   ! m
            model%climate%eus,                  &   ! m
            0.0d0,                              &   ! thklim = 0
            ice_mask,                           &
            floating_mask = floating_mask,      &
            ocean_mask = ocean_mask,            &
            land_mask = land_mask)

       ! Compute effective_areafrac for SMB purposes.
       ! Set = 1 where H > 0 so that thin ice will be removed from cells with a negative mass balance

       where (ice_mask == 1 .or. land_mask == 1)
          model%calving%effective_areafrac = 1.0d0
       elsewhere
          model%calving%effective_areafrac = 0.0d0
       endwhere

       if (verbose_calving) then
          call point_diag(model%calving%effective_areafrac, &
               'Before mass driver: effective_areafrac', itest, jtest, rtest, 7, 7, '(f10.6)')
       endif

    endif  ! which_ho_calving_front

    ! TODO: Zero out acab and bmlt in cells that are ice-free ocean after transport?
    !       Then it would not be necessary to pass ocean_mask to mass_balance_driver.

    ! Initialize the applied acab and bmlt.
    ! Note: These are smaller in magnitude than the input acab and bmlt for cells where either
    !       (1) the full column melts, and energy remains for melting, or
    !       (2) a positive mass balance is ignored, because a cell is ice-free ocean

    model%climate%acab_applied = 0.0d0
    model%basal_melt%bmlt_applied = 0.0d0

    ! ------------------------------------------------------------------------
    ! Apply the surface mass balance (acab) and basal mass balance (bmlt).
    ! Note: This subroutine assumes SI units:
    !       * dt (s)
    !       * dew, dns, thck (m)
    !       * acab, bmlt (m/s)
    ! ------------------------------------------------------------------------

    !TODO - Inline some of the code in this subroutine?
    !       As it is, this subroutine does very little other than call mass_balance_driver.

    call mass_balance_driver(&
         model%numerics%dt,                                    &  ! s
         model%numerics%dew,       model%numerics%dns,     &  ! m
         ewn,         nsn,         upn-1,                  &
         model%numerics%sigma,                             &
         parallel,                                         &
         itest,       jtest,       rtest,                  &
         model%geometry%thck(:,:),                         &  ! m
         model%climate%acab(:,:),                          &  ! m/s
         bmlt(:,:),                                        &  ! m/s
         model%climate%acab_applied(:,:),                  &  ! m/s
         model%basal_melt%bmlt_applied(:,:),               &  ! m/s
         ocean_mask(:,:),                                  &
         model%calving%effective_areafrac(:,:),            &
         model%geometry%ntracers,                          &
         model%geometry%tracers(:,:,:,:),                  &
         model%geometry%tracers_usrf(:,:,:),               &
         model%geometry%tracers_lsrf(:,:,:),               &
         model%options%which_ho_vertical_remap)

    !       End of mass balance code

  end subroutine glissade_apply_smb

!=======================================================================

  subroutine mass_balance_driver(&
       dt,                         &
       dx,           dy,           &
       nx,           ny,           &
       nlyr,         sigma,        &
       parallel,                   &
       itest, jtest, rtest,        &
       thck,                       &
       acab,         bmlt,         &
       acab_applied, bmlt_applied, &
       ocean_mask,                 &
       effective_areafrac,         &
       ntracers,     tracers,      &
       tracers_usrf, tracers_lsrf, &
       vert_remap_accuracy)

    ! This subroutine applies the surface and basal mass balance to each grid cell,
    !  keeping track of the total mass balance applied at each surface.
    ! It also checks for conservation of mass and mass*tracers.
    !
    ! Note: The SMB and BMB are not applied to ocean cells.
    !       For cells with an effective area less than 1 (e.g., calving-front cells),
    !        the SMB and BMB are applied only to the ice-covered part of the cell.
    !
    ! author William H. Lipscomb, NCAR
    !
    use glissade_transport, only: glissade_sum_mass_and_tracers, &
         glissade_vertical_remap, glissade_global_conservation

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

    real(dp), dimension(nx,ny), intent(in) ::  &
         effective_areafrac     ! effective fractional area, in range [0,1]
                                ! Calving_front cells can have values between 0 and 1

    real(dp), dimension(nx,ny), intent(inout) ::  &
         thck                   ! ice thickness (m), defined at horiz cell centers

    real(dp), dimension(nx,ny), intent(in) ::  &
         acab,    &             ! surface mass balance (m/s)
                                ! (defined at horiz cell centers)
         bmlt                   ! basal melt rate (m/s); positive for melting, negative for freeze-on
                                ! includes melting for both grounded and floating ice
                                ! (defined at horiz cell centers)

    ! Note: These fields are accumulated in units of meters, then converted to m/s.
    real(dp), dimension(nx,ny), intent(inout) ::  &
         acab_applied,    &     ! surface mass balance applied to ice (m)
                                ! = 0 for ice-free cells where acab < 0
         bmlt_applied           ! basal melt rate applied to ice (m)
                                ! = 0 for ice-free cells where bmlt > 0

    integer, dimension(nx,ny), intent(in) :: &
         ocean_mask             ! = 1 if topg is below sea level and thk <= thklim, else = 0

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

    ! local variables

    integer ::     &
         i, j, k         ,&! cell indices
         ilo,ihi,jlo,jhi ,&! beginning and end of physical domain
         nt                ! tracer index

    real(dp), dimension (nx,ny,nlyr) ::     &
         thck_layer        ! ice layer thickness (m)

    integer ::     &
         icells            ! number of cells with ice

    integer, dimension(nx*ny) ::     &
         indxi, indxj      ! compressed i/j indices

    real(dp) ::        &
         sum_acab,       & ! global sum of applied accumulation/ablation
         sum_bmlt          ! global sum of applied basal melting

    real(dp) ::     &
         msum_init,      &! initial global ice mass
         msum_final       ! final global ice mass

    !TODO - Delete these?
    real(dp), dimension(ntracers) ::     &
         mtsum_init,     &! initial global ice mass*tracer
         mtsum_final      ! final global ice mass*tracer

    real(dp), dimension(nx,ny) :: &
         melt_potential   ! total thickness (m) of additional ice that could be melted
                          ! by available acab/bmlt in columns that are completely melted

    real(dp) :: sum_melt_potential  ! global sum of melt potential

    logical ::     &
         errflag          ! true if energy is not conserved

    character(len=100) :: message

    real(dp) ::  &
         max_acab, max_bmlt  ! max magnitudes of acab and bmlt

    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

    errflag = .false.
    melt_potential(:,:) = 0.d0

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
    ! Add the mass balance at the surface and bed.
    ! Assume that new ice arrives at the surface with the current surface temperature.
    !
    ! TODO: Make sure this assumption is consistent with energy
    !       conservation for coupled simulations.
    ! TODO: Pass the melt potential back to the climate model as a heat flux?
    !-------------------------------------------------------------------

    max_acab = max(maxval(acab), -1.d0*minval(acab))
    max_bmlt = max(maxval(bmlt), -1.d0*minval(bmlt))

    max_acab = parallel_reduce_max(max_acab)
    max_bmlt = parallel_reduce_max(max_bmlt)

    if (max_acab > 0.0d0 .or. max_bmlt > 0.0d0) then

       call add_surface_and_basal_mass_balance(&
            nx,       ny,          &
            nlyr,     ntracers,    &
            dt,       parallel,    &
            ocean_mask,            &
            effective_areafrac,    &
            thck_layer(:,:,:),     &
            tracers(:,:,:,:),      &
            tracers_usrf(:,:,:),   &
            tracers_lsrf(:,:,:),   &
            acab(:,:),             &
            bmlt(:,:),             &
            acab_applied(:,:),     &
            bmlt_applied(:,:),     &
            melt_potential(:,:))

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
       ! Check that mass is conserved, allowing for mass gain/loss due to acab/bmlt
       !  and for any unused melt potential.
       !
       ! Note: There is no tracer conservation check here, because there is no
       !       easy way to correct initial mass*tracer values for acab and bmlt.
       !-------------------------------------------------------------------

       if (conservation_check) then

          ! Correct initial global mass for acab and bmlt
          sum_acab = 0.0d0
          sum_bmlt = 0.0d0
          sum_melt_potential = 0.0d0

          ! loop over locally owned cells, with correction for fractional coverage
          do j = 1+lhalo, ny-uhalo
             do i = 1+lhalo, nx-uhalo
                sum_acab = sum_acab + acab(i,j)*effective_areafrac(i,j)
                sum_bmlt = sum_bmlt + bmlt(i,j)*effective_areafrac(i,j)
                sum_melt_potential = sum_melt_potential + melt_potential(i,j)
             enddo
          enddo

          sum_acab = parallel_reduce_sum(sum_acab)
          sum_bmlt = parallel_reduce_sum(sum_bmlt)
          sum_melt_potential = parallel_reduce_sum(sum_melt_potential)

          msum_init = msum_init + (sum_acab - sum_bmlt)*dt

          ! Compute new global mass and mass*tracer

          call glissade_sum_mass_and_tracers(&
               nx,                ny,              &
               nlyr,              ntracers,        &
               thck_layer(:,:,:), msum_final,      &
               tracers(:,:,:,:),  mtsum_final(:))

          ! Check mass conservation
          !TODO - Add melt_potential to msum_final before calling subroutine?

          if (main_task) then

             call glissade_global_conservation(&
                  msum_init,     msum_final,      &
                  errflag,       sum_melt_potential)

             if (errflag) then
                write(message,*) 'WARNING: Conservation error in add_surface_and_basal_mass_balance'
!                  call write_log(message,GM_FATAL)      ! uncomment to make conservation errors fatal
                call write_log(message,GM_DIAGNOSTIC)  ! uncomment for debugging
             endif

          endif   ! main_task

       endif      ! conservation_check

       ! Recompute thickness
       thck(:,:) = 0.d0
       do k = 1, nlyr
          thck(:,:) = thck(:,:) + thck_layer(:,:,k)
       enddo

    endif  ! max_acab > 0 or max_bmlt > 0

  end subroutine mass_balance_driver

!=======================================================================

  subroutine add_surface_and_basal_mass_balance(&
       nx,           ny,          &
       nlyr,         ntracer,     &
       dt,           parallel,    &
       ocean_mask,                &
       effective_areafrac,        &
       thck_layer,   tracer,      &
       tracer_usrf,  tracer_lsrf, &
       acab,         bmlt,        &
       acab_applied, bmlt_applied,&
       melt_potential)

    ! Adjust the layer thickness based on the surface and basal mass balance

    ! Input/output arguments

    integer, intent(in) ::   &
         nx, ny,               &! horizontal array size
         nlyr,                 &! number of vertical layers
         ntracer                ! number of tracers

    real(dp), intent(in) ::   &
         dt                     ! time step (s)

    type(parallel_type), intent(in) :: &
         parallel               ! info for parallel communication

    !TODO - Could remove ocean_mask argument, if acab and bmlt have already been set to 0 for ice-free ocean cells.
    integer, dimension(nx,ny), intent(in) :: &
         ocean_mask             ! = 1 if topg is below sea level and thk <= thklim, else = 0

    real(dp), dimension (nx,ny), intent(in) ::     &
         effective_areafrac     ! effective fractional area for calving_front cells, in range [0,1]

    real(dp), dimension (nx,ny,nlyr), intent(inout) ::     &
         thck_layer             ! ice layer thickness

    real(dp), dimension (nx,ny,ntracer,nlyr), intent(inout) ::     &
         tracer                 ! 3D tracer values

    real(dp), dimension (nx,ny,ntracer), intent(in) ::     &
         tracer_usrf,          &! tracer values associated with accumulation at upper surface
         tracer_lsrf            ! tracer values associated with freeze-on at lower surface

    real(dp), intent(in), dimension(nx,ny) :: &
         acab                   ! surface mass balance (m/s)

    real(dp), intent(in), dimension(nx,ny) :: &
         bmlt                   ! basal melt rate (m/s)
                                ! > 0 for melting, < 0 for freeze-on

    real(dp), intent(inout), dimension(nx,ny) :: &
         acab_applied           ! surface mass balance applied to ice (m/s)
                                ! = 0 in ice-free regions where acab < 0

    real(dp), intent(inout), dimension(nx,ny) :: &
         bmlt_applied           ! basal melt rate applied to ice (m/s)
                                ! = 0 in ice-free regions where bmlt > 0

    real(dp), intent(out), dimension(nx,ny) :: &
         melt_potential   ! total thickness (m) of additional ice that could be melted
                          ! by available acab/bmlt in columns that are completely melted

    ! Local variables

    real(dp), dimension(nx,ny,ntracer,nlyr) ::  &
         thck_tracer       ! thck_layer * tracer

    real(dp), dimension(nx,ny) :: &
         thck_init,      & ! initial ice thickness
         thck_final        ! final ice thickness

    real(dp) :: sfc_accum, sfc_ablat  ! surface accumulation/ablation, from acab
    real(dp) :: bed_accum, bed_ablat  ! bed accumulation/ablation, from bmlt
    real(dp) :: dthck                 ! thickness change

    integer :: i, j, k, nt, iglobal, jglobal

    character(len=100) :: message

    ! Temporarily, convert the applied mass balance (intent inout) from m/s to m.
    ! It is converted back to m/s for output.
    acab_applied(:,:) = acab_applied(:,:) * dt
    bmlt_applied(:,:) = bmlt_applied(:,:) * dt

    ! Initialize the melt potential.
    ! These terms are adjusted below if energy is available for melting
    !  when no ice is present.

    melt_potential(:,:) = 0.0d0

    if (conservation_check) then
       do j = 1+nhalo, ny-nhalo
          do i = 1+nhalo, nx-nhalo
             thck_init(i,j) = sum(thck_layer(i,j,:))
          enddo
       enddo
    endif

    do j = 1+nhalo, ny-nhalo
       do i = 1+nhalo, nx-nhalo

          ! Temporarily adjust the layer thickness to account for partial ice converage.
          ! This prevents excessive thickening and thinning in partly filled calving front cells.
          if (effective_areafrac(i,j) > 0.0d0 .and. effective_areafrac(i,j) < 1.0d0) then
             thck_layer(i,j,:) = thck_layer(i,j,:) / effective_areafrac(i,j)
          endif

          ! initialize accumulation/ablation terms
          sfc_accum = 0.d0
          sfc_ablat = 0.d0
          bed_accum = 0.d0
          bed_ablat = 0.d0

          ! Add surface accumulation/ablation to ice thickness
          ! Also modify tracers conservatively.

          if (acab(i,j) > 0.d0) then       ! accumulation, added to layer 1

             sfc_accum = acab(i,j)*dt

             if (ocean_mask(i,j) == 1) then     ! no accumulation in open ocean

                ! do nothing

             else  ! not ocean; accumulate ice

                acab_applied(i,j) = acab_applied(i,j) + sfc_accum*effective_areafrac(i,j)

                ! adjust mass-tracer product for the top layer

                do nt = 1, ntracer  !TODO - Put this loop on the outside for speedup?

                   thck_tracer(i,j,nt,1) = thck_layer(i,j,1) * tracer(i,j,nt,1)  &
                                                 + sfc_accum * tracer_usrf(i,j,nt)

                enddo  ! ntracer

                ! new top layer thickess
                thck_layer(i,j,1) = thck_layer(i,j,1) + sfc_accum

                ! new tracer values in top layer
                tracer(i,j,:,1) = thck_tracer(i,j,:,1) / thck_layer(i,j,1)

             endif   ! ocean_mask = 1

          elseif (acab(i,j) < 0.d0) then   ! ablation in one or more layers

             ! reduce ice thickness (tracer values will not change)

             sfc_ablat = -acab(i,j)*dt   ! positive by definition

             if (ocean_mask(i,j) == 1) then     ! no accumulation in open ocean

                ! do nothing

             else  ! not ocean; melt ice

                acab_applied(i,j) = acab_applied(i,j) - sfc_ablat*effective_areafrac(i,j)

                do k = 1, nlyr
                   if (sfc_ablat > thck_layer(i,j,k)) then
                      sfc_ablat = sfc_ablat - thck_layer(i,j,k)
                      thck_layer(i,j,k) = 0.d0
                      tracer(i,j,:,k) = 0.d0
                   else
                      thck_layer(i,j,k) = thck_layer(i,j,k) - sfc_ablat
                      sfc_ablat = 0.d0
                      exit
                   endif
                enddo

             endif   ! ocean_mask = 1

             ! Adjust acab_applied if energy is still available for melting
             ! Also accumulate the remaining melt energy

             if (sfc_ablat > 0.d0) then
                acab_applied(i,j) = acab_applied(i,j) + sfc_ablat*effective_areafrac(i,j)  ! make a negative value less negative
                melt_potential(i,j) = melt_potential(i,j) + sfc_ablat
             endif

             !TODO - Figure out how to handle excess energy given by melt_potential.
             !       Include in the heat flux passed back to CLM?

          endif  ! acab > 0

          ! Note: It is possible that we could have residual energy remaining for surface ablation
          !       while ice is freezing on at the bed, in which case the surface ablation should
          !       be subtracted from the bed accumulation.  We ignore this possibility for now.

          ! Note: Freeze-on (bmlt < 0) is allowed only in ice-covered cells, not ice-free ocean.
          !       Allowing freeze-on in ice-free ocean would introduce mass conservation errors,
          !        given the current logic with effective_areafrac.
          !       If it is desired to implement a field of frazil ice formation that could grow ice
          !        in open ocean as well as sub-shelf cavities and open ocean, this field could be passed
          !        into glissade_mass_balance_driver in a separate call (i.e., independent of the standard
          !        acab and bmlt fields) with ocean_mask = 0 and effective_areafrac = 1 everywhere.
          !        Then the following code would allow frazil growth, as desired.

          if (bmlt(i,j) < 0.d0) then       ! freeze-on, added to lowest layer

             bed_accum = -bmlt(i,j)*dt

             if (ocean_mask(i,j) == 1) then     ! no accumulation in open ocean

                ! do nothing

             else  ! not ocean; accumulate ice

                bmlt_applied(i,j) = bmlt_applied(i,j) - bed_accum*effective_areafrac(i,j)  ! bmlt_applied < 0 for freeze-on

                ! adjust mass-tracer product for the bottom layer

                do nt = 1, ntracer  !TODO - Put this loop on the outside for speedup?

                   thck_tracer(i,j,nt,nlyr) = thck_layer(i,j,nlyr) * tracer(i,j,nt,nlyr)  &
                                                       + bed_accum * tracer_lsrf(i,j,nt)

                enddo  ! ntracer

                ! new bottom layer thickess
                thck_layer(i,j,nlyr) = thck_layer(i,j,nlyr) + bed_accum

                ! new tracer values in bottom layer
                tracer(i,j,:,nlyr) = thck_tracer(i,j,:,nlyr) / thck_layer(i,j,nlyr)

             endif   ! ocean_mask = 1

          elseif (bmlt(i,j) > 0.d0) then   ! basal melting in one or more layers

             ! reduce ice thickness (tracer values will not change)

             bed_ablat = bmlt(i,j)*dt   ! positive by definition

             if (ocean_mask(i,j) == 1) then     ! no accumulation in open ocean

                ! do nothing

             else  ! not ocean; melt ice

                bmlt_applied(i,j) = bmlt_applied(i,j) + bed_ablat*effective_areafrac(i,j)

                do k = nlyr, 1, -1
                   if (bed_ablat > thck_layer(i,j,k)) then
                      bed_ablat = bed_ablat - thck_layer(i,j,k)
                      thck_layer(i,j,k) = 0.d0
                      tracer(i,j,:,k) = 0.d0
                   else
                      thck_layer(i,j,k) = thck_layer(i,j,k) - bed_ablat
                      bed_ablat = 0.d0
                      exit
                   endif
                enddo

             endif   ! ocean_mask = 1

             ! Adjust bmlt_applied if energy is still available for melting
             ! Also accumulate the remaining melt energy

             if (bed_ablat > 0.d0) then
                ! bmlt_applied is less than input bmlt
                bmlt_applied(i,j) = bmlt_applied(i,j) - bed_ablat*effective_areafrac(i,j)
                melt_potential(i,j) = melt_potential(i,j) + bed_ablat
             endif

          endif  ! bmlt < 0

          ! Weight the melt potential by the effective area fraction
          melt_potential(i,j) = melt_potential(i,j) * effective_areafrac(i,j)

          ! Convert thck_layer back to the mean volume per unit area in partly covered cells
          if (effective_areafrac(i,j) > 0.0d0 .and. effective_areafrac(i,j) < 1.0d0) then
             thck_layer(i,j,:) = thck_layer(i,j,:) * effective_areafrac(i,j)
          endif

       enddo   ! i
    enddo      ! j

    ! Check mass conservation in each column

    if (conservation_check) then
       do j = 1+nhalo, ny-nhalo
          do i = 1+nhalo, nx-nhalo
             thck_final(i,j) = sum(thck_layer(i,j,:))
             dthck = (acab(i,j) - bmlt(i,j))*dt*effective_areafrac(i,j)
             if (abs(thck_init(i,j) + dthck - thck_final(i,j) + melt_potential(i,j)) > 1.d-8) then
                call parallel_globalindex(i, j, iglobal, jglobal, parallel)
                write(6,*) ' '
                write(6,*) 'ERROR: Column conservation check, r, i, j, iglobal, jglobal, err =', &
                     this_rank, i, j, iglobal, jglobal, thck_init(i,j) + dthck - thck_final(i,j)
                write(6,*) 'thck_init, dthck, thck_final:', thck_init(i,j), dthck, thck_final(i,j)
                write(6,*) 'acab*dt, bmlt*dt, areafrac, melt_potential:', &
                     acab(i,j)*dt, bmlt(i,j)*dt, effective_areafrac(i,j), melt_potential(i,j)
                write(message,*) &
                     'WARNING: Column conservation error in add_surface_and_basal_mass_balance, i, j =', i, j
!                call write_log(message,GM_FATAL)       ! uncomment to make conservation errors fatal
                call write_log(message,GM_DIAGNOSTIC)  ! uncomment for debugging
             endif
          enddo
       enddo
    endif

    ! convert applied mass balance from m to m/s
    acab_applied(:,:) = acab_applied(:,:) / dt
    bmlt_applied(:,:) = bmlt_applied(:,:) / dt

  end subroutine add_surface_and_basal_mass_balance

!=======================================================================

  subroutine set_overwrite_acab_mask(&
       overwrite_acab,         &
       acab,                   &
       thck,                   &
       overwrite_acab_minthck, &
       overwrite_acab_mask)

    use glide_types

    ! If overwrite_acab /=0 , then set overwrite_acab_mask = 1 for grid cells
    !  where acab is to be overwritten.  Currently, three options are supported:
    ! (1) Overwrite acab where the input acab = 0 at initialization
    ! (2) Overwrite acab where the input thck <= overwrite_acab_minthck at initialization
    ! (3) Overwrite acab based on an input mask
    !
    ! Note: This subroutine should be called only on initialization, not on restart.

    integer, intent(in) ::  &
         overwrite_acab           !> option for overwriting acab

    real(dp), dimension(:,:), intent(in) ::  &
         acab,                  & !> ice surface mass balance (model units)
         thck                     !> ice thickness (model units)

    real(dp), intent(in) ::  &
         overwrite_acab_minthck   !> overwrite acab where thck <= overwrite_acab_minthck (model units)

    integer, dimension(:,:), intent(out) ::  &
         overwrite_acab_mask      !> = 1 where acab is overwritten, else = 0

    integer :: ewn, nsn
    integer :: i, j
    integer :: max_mask_local, max_mask_global

    ewn = size(overwrite_acab_mask,1)
    nsn = size(overwrite_acab_mask,2)

    if (overwrite_acab == OVERWRITE_ACAB_ZERO_ACAB) then

       do j = 1, nsn
          do i = 1, ewn

             if (acab(i,j) == 0.0d0) then
                overwrite_acab_mask(i,j) = 1
             else
                overwrite_acab_mask(i,j) = 0
             endif

          enddo
       enddo

    elseif (overwrite_acab == OVERWRITE_ACAB_THCKMIN) then

       do j = 1, nsn
          do i = 1, ewn

             ! Note the '<='.  If overwrite_acab_minthck = 0.d0, only ice-free cells are overwritten.
             if (thck(i,j) <= overwrite_acab_minthck) then
                overwrite_acab_mask(i,j) = 1
             else
                overwrite_acab_mask(i,j) = 0
             endif

          enddo
       enddo

    elseif (overwrite_acab == OVERWRITE_ACAB_INPUT_MASK) then

       ! Make sure a mask was read in with some nonzero values
       ! If not, then write a warning

       max_mask_local = maxval(overwrite_acab_mask)
       max_mask_global = parallel_reduce_max(max_mask_local)
       if (main_task) then
          write(6,*) 'rank, max_mask_local, max_mask_global:', &
               this_rank, max_mask_local, max_mask_global
       endif
       if (max_mask_global == 1) then
          ! continue
       elseif (max_mask_global == 0) then
          call write_log('Using overwrite_acab_mask without any values > 0', GM_WARNING)
       else
          call write_log('Using overwrite_acab_mask with values other than 0 and 1', GM_FATAL)
       endif

    endif  ! overwrite_acab

  end subroutine set_overwrite_acab_mask

!=======================================================================

  subroutine overwrite_acab(&
       overwrite_acab_mask,  &
       overwrite_acab_value, &
       acab)

    integer, dimension(:,:), intent(in) ::  &
         overwrite_acab_mask     !> mask = 1 where acab is overwritten value, else = 0

    real(dp), intent(in) ::  &
         overwrite_acab_value    !> acab value applied where overwrite_acab_mask = 1

    real(dp), dimension(:,:), intent(inout) ::  &
         acab           !> unadjusted acab (model units) on input
                        !> overwritten acab on output

    integer :: ewn, nsn
    integer :: i, j

    ewn = size(acab,1)
    nsn = size(acab,2)

    do j = 1, nsn
       do i = 1, ewn

          if (overwrite_acab_mask(i,j) == 1) then
             acab(i,j) = overwrite_acab_value
          endif

       enddo
    enddo

  end subroutine overwrite_acab

!=======================================================================

  subroutine glissade_add_2d_anomaly(&
       var2d,                    &
       var2d_anomaly,            &
       anomaly_tstart,           &
       anomaly_timescale,        &
       time)

    ! Apply a 2D anomaly field, usually to the surface mass balance, surface temperature,
    ! or similar climate forcing field.

    use glimmer_paramets, only: eps08

    real(dp), dimension(:,:), intent(inout) ::  &
         var2d               !> 2D field (uncorrected)
                             !> uncorrrected on input, corrected on output

    real(dp), dimension(:,:), intent(in) ::   &
         var2d_anomaly       !> anomalous field to be added to the var2d input value

    real(dp), intent(in) ::  &
         anomaly_tstart,   & !> time to begin applying the anomaly (yr)
         anomaly_timescale   !> number of years over which the anomaly is phased in linearly

    real(dp), intent(in) :: &
         time                !> model time in years
                             !> Note: Should be the time at the start of the time step, not the end

    integer :: ewn, nsn
    integer :: i, j
    real(dp) :: anomaly_fraction

    ewn = size(var2d,1)
    nsn = size(var2d,2)

    ! Given the model time, compute the fraction of the anomaly to be applied now.
    ! Add a small value to the time to avoid rounding errors when time is close to an integer value.

    if (time + eps08 > anomaly_tstart + anomaly_timescale .or. anomaly_timescale == 0.0d0) then

       ! apply the full anomaly
       anomaly_fraction = 1.0d0

    elseif (time + eps08 > anomaly_tstart) then

       ! apply an increasing fraction of the anomaly
       ! Note: There are three options:
       ! (1) Apply the anomaly in proportion to the elapsed time since anomaly_start
       ! (2) Apply the anomaly in one-year chunks, starting one year after anomaly_start.
       !     This is the convention for initMIP.
       ! (3) Apply the anomaly in one-year chunks, starting immediately after anomaly_start.
       !     That is, the one-year anomaly is applied throughout the first year.
       ! For now, option (3) is the default.

!       anomaly_fraction = (time - anomaly_tstart) / anomaly_timescale
!       anomaly_fraction = floor(time + eps08 - anomaly_tstart, dp) / anomaly_timescale

       ! Subtract a small number from the time, given that model%numerics%time is the time
       ! at the end of the timestep, so e.g. the December time at the end of 1983 would be
       ! close to 1984.0, and we don't want to jump to the following year due to rounding error.

       anomaly_fraction = ceiling(time - eps08 - anomaly_tstart, dp) / anomaly_timescale
!!       if (main_task) write(6,*) 'In add_2d_anomaly: time, frac =', time, anomaly_fraction         

    else
       ! no anomaly to apply
       anomaly_fraction = 0.0d0
    endif

    ! apply the anomaly
    do j = 1, nsn
       do i = 1, ewn
          var2d(i,j) = var2d(i,j) + anomaly_fraction*var2d_anomaly(i,j)
       enddo
    enddo

  end subroutine glissade_add_2d_anomaly

!=======================================================================

  !TODO - Make the logic consistent with the 2d subroutine above
  subroutine glissade_add_3d_anomaly(var3d,                 &
                                     var3d_anomaly,         &
                                     anomaly_tstart,        &
                                     anomaly_timescale,     &
                                     time)

    ! Apply a 3D anomaly field, usually to the surface mass balance, surface temperature,
    ! or similar climate forcing field. Not currently called.

    use glimmer_paramets, only : eps08

    real(dp), dimension(:,:,:), intent(inout) ::  &
         var3d               !> input 3d field; vertical index in first slot
                             !> uncorrrected on input; anomaly added on output

    real(dp), dimension(:,:,:), intent(in) ::   &
         var3d_anomaly       !> anomaly to be added to the input value

    real(dp), intent(in) ::  &
         anomaly_tstart,   & !> time to begin applying the anomaly (yr)
         anomaly_timescale   !> number of years over which the anomaly is phased in linearly

    real(dp), intent(in) :: &
         time                     !> model time in years
                                  !> Note: Should be the time at the start of the time step, not the end

    integer :: ewn, nsn
    integer :: i, j
    real(dp) :: anomaly_fraction

    ewn = size(var3d,2)
    nsn = size(var3d,3)

    ! Given the model time, compute the fraction of the anomaly to be applied now.
    ! Add a small value to the time to avoid rounding errors when time is close to an integer value.

    if (time + eps08 > anomaly_tstart + anomaly_timescale .or. anomaly_timescale == 0.0d0) then

       ! apply the full anomaly
       anomaly_fraction = 1.0d0

    elseif (time + eps08 > anomaly_tstart) then

       ! apply an increasing fraction of the anomaly
       anomaly_fraction = (time - anomaly_tstart) / anomaly_timescale

       ! Note: For initMIP, the anomaly is applied in annual step functions
       !        starting at the end of the first year.
       !       Comment out the line above and uncomment the following line
       !        to increase the anomaly once a year.
!       anomaly_fraction = floor(time + eps08 - anomaly_tstart, dp) / anomaly_timescale
!
    else
       ! no anomaly to apply
       anomaly_fraction = 0.0d0
    endif

    ! apply the anomaly

    do j = 1, nsn
       do i = 1, ewn
          var3d(:,i,j) = var3d(:,i,j) + anomaly_fraction*var3d_anomaly(:,i,j)
       enddo
    enddo

  end subroutine glissade_add_3d_anomaly

!=======================================================================

  end module glissade_mass_balance

!=======================================================================
