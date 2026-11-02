!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_diagnostics.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glide_diagnostics

  ! subroutines for computing various useful diagnostics
  ! Author: William Lipscomb

  use glimmer_global, only: dp
  use glimmer_paramets, only: iulog, eps11
  use glimmer_log
  use glide_types

  ! Note: Glide runs on a single task with no halo cells, so these diagnostics
  !       do not need parallel reductions or broadcasts.

  implicit none

  private :: masked_sum_2d, masked_sum_3d

  logical :: verbose_diagnostics = .false.

contains

  subroutine glide_write_diagnostics (model, time, tstep_count)

    ! Short driver subroutine to decide whether it's time to write diagnostics.
    ! If so, it calls glide_write_diag.

    ! input/output arguments

    ! Note: model is intent(inout) so that some global diagnostics can be computed below
    type(glide_global_type), intent(inout) :: model    ! model instance

    real(dp), intent(in) :: time          ! current time in years

    integer,  intent(in) :: tstep_count   ! current timestep

    ! local arguments

    ! debug
    if (verbose_diagnostics) then
       write(iulog,*) '	'
       write(iulog,*) 'In glide_write_diagnostics'
       write(iulog,*) 'time =', time
       write(iulog,*) 'dt_diag =', model%numerics%dt_diag
       write(iulog,*) 'ndiag =', model%numerics%ndiag
       write(iulog,*) 'tstep_count =', tstep_count
    endif

    if (model%numerics%ndiag > 0) then

       if (mod(tstep_count, model%numerics%ndiag) == 0)  then    ! time to write
          call glide_write_diag(model, time)
       endif

    endif    ! ndiag > 0

  end subroutine glide_write_diagnostics
 
!--------------------------------------------------------------------------

  subroutine glide_init_diag (model)

    ! Initialize model diagnostics for Glide.
    ! (1) Set ndiag based on dt_diag.  (Diagnostics are written every ndiag steps.)
    ! (2) Set the local indices and rank of the global diagnostic point


    implicit none

    ! input/output arguments

    type(glide_global_type), intent(inout) :: model    ! model instance

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
    ! Set the local indices and rank of the global diagnostic point.
    ! Glide runs on a single task (rank 0) with no halo cells, so the local
    ! indices are the same as the global indices.
    ! Note: The global indices are checked in glide_setup (print_parameters).
    !-----------------------------------------------------------------

    model%numerics%idiag_local = model%numerics%idiag
    model%numerics%jdiag_local = model%numerics%jdiag
    model%numerics%rdiag_local = 0

    write(iulog,'(a26,2i6)') ' Global idiag, jdiag:     ',   &
                          model%numerics%idiag, model%numerics%jdiag
    write(iulog,'(a26,3i6)') ' Local idiag, jdiag, task:',   &
                          model%numerics%idiag_local,  &
                          model%numerics%jdiag_local,  &
                          model%numerics%rdiag_local

    write(message,'(a25,2i6)') 'Global idiag, jdiag:     ',   &
                                model%numerics%idiag, model%numerics%jdiag
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,3i6)') 'Local idiag, jdiag, task:',   &
                                model%numerics%idiag_local,  &
                                model%numerics%jdiag_local,  &
                                model%numerics%rdiag_local
    call write_log(trim(message), type = GM_DIAGNOSTIC)

  end subroutine glide_init_diag

!--------------------------------------------------------------------------

  subroutine glide_write_diag (model, time)

    ! Write global diagnostics
    ! Also write local diagnostics for a selected grid cell
 
    use glimmer_paramets, only: unphys_val
    use glimmer_physcon, only: scyr, rhoi, shci

    implicit none
 
    ! input/output arguments

    ! Note: model is intent(inout) so that some global diagnostics can be computed here 
    type(glide_global_type), intent(inout) :: model ! model instance

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
         thck_floating,                 &    ! thickness of floating ice
         thck_above_flotation,          &    ! thickness above flotation
         tot_energy,                    &    ! total ice energy (J)
         mean_thck,                     &    ! mean ice thickness (m)
         mean_temp,                     &    ! mean ice temperature (deg C)
         max_thck,                      &    ! max ice thickness (m)
         max_temp,                      &    ! max ice temperature (deg C)
         min_temp,                      &    ! min ice temperature (deg C)
         max_bmlt,                      &    ! max basal melt rate (m/yr)
         max_spd_sfc,                   &    ! max surface ice speed (m/yr)
         max_spd_bas,                   &    ! max basal ice speed (m/yr)
         spd,                           &    ! speed
         thck_diag, usrf_diag,          &    ! local column diagnostics
         topg_diag, relx_diag,          &    
         load_diag,                     &
         artm_diag, acab_diag,          &
         bmlt_diag, bwat_diag,          &
         bheatflx_diag,                 &
         level

    integer, dimension(model%general%ewn,model%general%nsn) ::  &
         ice_mask,                 & ! = 1 where ice is present with thck > minthck, else = 0
         floating_mask,            & ! = 1 where ice is present and floating, else = 0
         grounded_mask               ! = 1 where ice is present and grounded, else = 0

    integer, dimension(model%general%ewn-1,model%general%nsn-1) ::  &
         stag_ice_mask               ! staggered mask; = 1 if ice_mask = 1 for any of the four neighbors

    real(dp), dimension(model%general%upn) ::  &
         temp_diag,                     &    ! Note: sfc temp not included if temps are staggered
                                             !       (use artm instead)
         spd_diag                            ! speed (m/yr)

    real(dp), dimension(model%lithot%nlayer) ::  &
         lithtemp_diag                       ! lithosphere column diagnostics

    integer :: &
         i, j, k,                       &
         ktop, kbed,                        &
         imax, imin,                        &
         jmax, jmin,                        &
         kmax, kmin,                        &
         ewn, nsn, upn,                     &    ! model%numerics%ewn, etc.
         nlith                                   ! model%lithot%nlayer

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
         volume_above_flotation   ! ice volue above flotation (m^3)

    real(dp), dimension(model%general%upn, model%general%ewn, model%general%nsn) ::  &
         local_energy           ! internal energy (J) per layer in a column

    character(len=100) :: message
    
    ! Note: cell_area is copied here from model%geometry%cell_area
    ! cell_area = dew*dns by default; optionally scaled to account for grid distortion
    real(dp), dimension(model%general%ewn,model%general%nsn) :: &
         cell_area     ! grid cell areas (scaled model units); diagnostic only

    real(dp), parameter ::   &
         eps = eps11             ! small threshold for diagnostics

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

    !-----------------------------------------------------------------
    ! Compute and write global diagnostics
    !-----------------------------------------------------------------
 
    call write_log('----------------------------------------------------------')
    call write_log(' ')
    write(message,'(a25,f24.16)') 'Diagnostic output, time =', time
    call write_log(trim(message), type = GM_DIAGNOSTIC)
    call write_log(' ')

    ! total ice area (m^2)
    tot_area = masked_sum_2d(cell_area, ice_mask)
    tot_area_float = masked_sum_2d(cell_area, floating_mask)
    tot_area_ground = masked_sum_2d(cell_area, grounded_mask)

    ! total ice volume (m^3)
    tot_volume = masked_sum_2d(model%geometry%thck*cell_area, ice_mask)

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

    tot_volume_above_flotation = masked_sum_2d(volume_above_flotation)

    ! total ice mass above flotation (kg)
    tot_mass_above_flotation = tot_volume_above_flotation * rhoi

    ! total ice energy relative to T = 0 deg C (J)
    local_energy = 0.0d0
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
    local_energy = local_energy * rhoi * shci
    tot_energy = masked_sum_3d(local_energy, ice_mask)

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

    write(message,'(a25,f24.16)') 'Mean thickness (m)       ', mean_thck
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,f24.16)') 'Mean temperature (C)     ', mean_temp
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! Find various global maxes and mins

    ! max thickness

    imax = 0
    jmax = 0
    max_thck = unphys_val   ! = an arbitrary large negative number
    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > max_thck) then
             max_thck = model%geometry%thck(i,j)
             imax = i
             jmax = j
          endif
       enddo
    enddo

    write(message,'(a25,f24.16,2i6)') 'Max thickness (m), i, j  ',   &
                                       max_thck, imax, jmax
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max temperature

    ktop = lbound(model%temper%temp,1)
    kbed = ubound(model%temper%temp,1)

    imax = 0
    jmax = 0
    kmax = 0
    max_temp =  unphys_val
    do j = 1, nsn
       do i = 1, ewn
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

    write(message,'(a25,f24.16,3i6)') 'Max temperature, i, j, k ',   &
                    max_temp, imax, jmax, kmax
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! min temperature

    imin = 0
    jmin = 0
    kmin = 0
    min_temp = 999.d0  ! arbitrary large positive number
    do j = 1, nsn
       do i = 1, ewn
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

    write(message,'(a25,f24.16,3i6)') 'Min temperature, i, j, k ',   &
                    min_temp, imin, jmin, kmin
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max applied basal melt rate
    ! Usually, this will be for floating ice, if floating ice is present
    imax = 0
    jmax = 0
    max_bmlt = unphys_val

    do j = 1, nsn
       do i = 1, ewn
          if (model%basal_melt%bmlt_applied(i,j) > max_bmlt) then
             max_bmlt = model%basal_melt%bmlt_applied(i,j)
             imax = i
             jmax = j
          endif
       enddo
    enddo

    ! Write to diagnostics only if nonzero

    if (abs(max_bmlt*scyr) > eps) then
       write(message,'(a25,f24.16,2i6)') 'Max bmlt (m/y), i, j     ',   &
            max_bmlt*scyr, imax, jmax
       call write_log(trim(message), type = GM_DIAGNOSTIC)
    endif

    ! max surface speed
    imax = 0
    jmax = 0
    max_spd_sfc = unphys_val

    do j = 1, nsn-1          ! velocity grid
       do i = 1, ewn-1
          spd = sqrt(model%velocity%uvel(1,i,j)**2   &
                   + model%velocity%vvel(1,i,j)**2)
          if (model%geometry%stagthck(i,j) > minthck .and. spd > max_spd_sfc) then
             max_spd_sfc = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo

    write(message,'(a25,f24.16,2i6)') 'Max sfc spd (m/y), i, j  ',   &
                    max_spd_sfc*scyr, imax, jmax
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max basal speed

    imax = 0
    jmax = 0
    max_spd_bas = unphys_val
    do j = 1, nsn-1          ! velocity grid
       do i = 1, ewn-1
          spd = sqrt(model%velocity%uvel(upn,i,j)**2   &
                   + model%velocity%vvel(upn,i,j)**2)
          if (model%geometry%stagthck(i,j) > minthck  .and. spd > max_spd_bas) then
             max_spd_bas = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo

    write(message,'(a25,f24.16,2i6)') 'Max base spd (m/y), i, j ',   &
                    max_spd_bas*scyr, imax, jmax
    call write_log(trim(message), type = GM_DIAGNOSTIC)

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
    temp_diag(:)  = unphys_val
    spd_diag(:)   = unphys_val
    lithtemp_diag(:) = unphys_val    

    ! Set local diagnostic values

    if (model%numerics%idiag_local >= 1 .and. model%numerics%idiag_local <= ewn  &
                                        .and.                                    &
        model%numerics%jdiag_local >= 1 .and. model%numerics%jdiag_local <= nsn) then

       i = model%numerics%idiag_local
       j = model%numerics%jdiag_local
       usrf_diag = model%geometry%usrf(i,j)
       thck_diag = model%geometry%thck(i,j)
       topg_diag = model%geometry%topg(i,j)
       if (model%options%isostasy == ISOSTASY_COMPUTE) then
          relx_diag = model%isostasy%relx(i,j)
          load_diag = model%isostasy%load(i,j)
       endif
       artm_diag = model%climate%artm(i,j)
       acab_diag = model%climate%acab_applied(i,j) * scyr
       bmlt_diag = model%basal_melt%bmlt_applied(i,j) * scyr
       bwat_diag = model%basal_hydro%bwat(i,j)
       bheatflx_diag = model%temper%bheatflx(i,j)
       temp_diag(:) = model%temper%temp(1:upn,i,j)
       spd_diag(:) = sqrt(model%velocity%uvel(1:upn,i,j)**2   &
                        + model%velocity%vvel(1:upn,i,j)**2) * scyr
       if (model%options%gthf == GTHF_COMPUTE) &
            lithtemp_diag(:) = model%lithot%temp(i,j,:)

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

       ! Vertical profile of ice speed and temperature

       call write_log(' ')
       write(message,'(a)') ' Sigma       Ice speed (m/yr)      Ice temperature (C)'
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       ! temperatures are at layer interfaces, with speeds
       do k = 1, upn
          write (message,'(f6.3,2f24.16)') model%numerics%sigma(k), spd_diag(k), temp_diag(k)
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       enddo

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

  end subroutine glide_write_diag

!==============================================================

  function masked_sum_2d(a, mask)

    ! Sum of a 2D field over cells with mask = 1 (all cells if mask is absent).
    ! Note: The loop order matches parallel_global_sum (without reprosum),
    !       which these diagnostics used before.

    real(dp), dimension(:,:), intent(in) :: a
    integer, dimension(:,:), intent(in), optional :: mask
    real(dp) :: masked_sum_2d

    integer :: i, j

    masked_sum_2d = 0.0d0
    do j = 1, size(a,2)
       do i = 1, size(a,1)
          if (present(mask)) then
             if (mask(i,j) /= 1) cycle
          endif
          masked_sum_2d = masked_sum_2d + a(i,j)
       enddo
    enddo

  end function masked_sum_2d

!==============================================================

  function masked_sum_3d(a, mask)

    ! Sum of a 3D field a(k,i,j) over columns with mask(i,j) = 1.
    ! Note: The loop order matches parallel_global_sum (without reprosum),
    !       which these diagnostics used before.

    real(dp), dimension(:,:,:), intent(in) :: a
    integer, dimension(:,:), intent(in) :: mask
    real(dp) :: masked_sum_3d

    integer :: i, j, k

    masked_sum_3d = 0.0d0
    do j = 1, size(a,3)
       do i = 1, size(a,2)
          if (mask(i,j) == 1) then
             do k = 1, size(a,1)
                masked_sum_3d = masked_sum_3d + a(k,i,j)
             enddo
          endif
       enddo
    enddo

  end function masked_sum_3d

!==============================================================

end module glide_diagnostics

!==============================================================
