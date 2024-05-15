!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WARNING: this file was automatically generated on
! Tue, 07 May 2024 09:12:04 +0000
! from ncdf_template.F90.in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! WJS (1-30-12): The following (turning optimization off) is needed as a workaround for an
! xlf compiler bug, at least in IBM XL Fortran for AIX, V12.1 on bluefire
#ifdef CPRIBM
@PROCESS OPT(0)
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   ncdf_template.F90.in - part of the Community Ice Sheet Model (CISM)  
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

#define NCO outfile%nc
#define NCI infile%nc

#define HAVE_AVG 1

module glide_io
  ! template for creating subsystem specific I/O routines
  ! written by Magnus Hagdorn, 2004

  !WHL, Sept. 2023
  ! Moved some 'use' statements to the top to avoid redundant statements that slow the build,
  ! particularly on the Intel compiler

  use glide_types
  use glimmer_ncdf
  use glimmer_paramets
  use glimmer_physcon
  use glimmer_scales

  implicit none

  private :: get_xtype, is_enabled, is_enabled_0dint, is_enabled_1dint, &
       is_enabled_2dint, is_enabled_0dreal, is_enabled_1dreal, is_enabled_2dreal, is_enabled_3dreal

  character(glimmer_nc_vars_len), save :: restart_variable_list=''    ! list of variables needed for a restart

  interface is_enabled  ! MJH 10/21/13: Interface needed for determining if arrays have been enabled.  See notes below in glide_io_create.
    module procedure is_enabled_0dint
    module procedure is_enabled_1dint
    module procedure is_enabled_2dint
    module procedure is_enabled_0dreal
    module procedure is_enabled_1dreal
    module procedure is_enabled_2dreal
    module procedure is_enabled_3dreal
  end interface is_enabled

contains

  !*****************************************************************************
  ! netCDF output
  !*****************************************************************************
  subroutine glide_io_createall(model,data,outfiles)
    ! open all netCDF files for output
    use glimmer_ncio
    implicit none
    type(glide_global_type) :: model
    type(glide_global_type) :: data ! MJH 10/21/13: Making 'data' mandatory.  See notes below in glide_io_create
    type(glimmer_nc_output),optional,pointer :: outfiles
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    do while(associated(oc))
       call glide_io_create(oc,model,data)
       oc=>oc%next
    end do
  end subroutine glide_io_createall

  subroutine glide_io_writeall(data,model,atend,outfiles,time)
    ! if necessary write to netCDF files
    use glimmer_ncio
    implicit none
    type(glide_global_type) :: data
    type(glide_global_type) :: model
    logical, optional :: atend
    type(glimmer_nc_output),optional,pointer :: outfiles
    real(dp),optional :: time

    ! local variables
    type(glimmer_nc_output), pointer :: oc
    logical :: forcewrite=.false.

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    if (present(atend)) then
       forcewrite = atend
    end if

    do while(associated(oc))
#ifdef HAVE_AVG
       if (oc%do_averages) then
          call glide_avg_accumulate(oc,data,model)
       end if
#endif
       call glimmer_nc_checkwrite(oc,model,forcewrite,time)
       if (oc%nc%just_processed) then
          ! write standard variables
          call glide_io_write(oc,data)
#ifdef HAVE_AVG
          if (oc%do_averages) then
             call glide_avg_reset(oc,data)
          end if
#endif
       end if
       oc=>oc%next
    end do
  end subroutine glide_io_writeall
  
  subroutine glide_io_create(outfile,model,data)

    use cism_parallel, only: parallel_type, &
         parallel_def_dim, parallel_inq_dimid, parallel_def_var, parallel_inq_varid, parallel_put_att
    use glimmer_ncio
    use glimmer_map_types
    use glimmer_log
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    type(glide_global_type) :: model
    type(glide_global_type) :: data    ! MJH 10/21/13: Making 'data' mandatory.  See note below

    integer status,varid,pos

    ! MJH 10/21/13: Local variables needed for checking if a variable is enabled.
    real(dp) :: tavgf
    integer :: up

    integer :: level_dimid
    integer :: lithoz_dimid
    integer :: nlev_smb_dimid
    integer :: staglevel_dimid
    integer :: stagwbndlevel_dimid
    integer :: time_dimid
    integer :: x0_dimid
    integer :: x1_dimid
    integer :: y0_dimid
    integer :: y1_dimid
    integer :: zocn_dimid

    ! defining dimensions
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'level',model%general%upn,level_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'level',level_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'lithoz',model%lithot%nlayer,lithoz_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'lithoz',lithoz_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'nlev_smb',model%climate%nlev_smb,nlev_smb_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'nlev_smb',nlev_smb_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'staglevel',model%general%upn-1,staglevel_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'staglevel',staglevel_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'stagwbndlevel',model%general%upn+1,stagwbndlevel_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'stagwbndlevel',stagwbndlevel_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inq_dimid(NCO%id,'time',time_dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'x0',model%parallel%global_ewn-1,x0_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'x0',x0_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'x1',model%parallel%global_ewn,x1_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'x1',x1_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'y0',model%parallel%global_nsn-1,y0_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'y0',y0_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'y1',model%parallel%global_nsn,y1_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'y1',y1_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (.not.outfile%append) then
       status = parallel_def_dim(NCO%id,'zocn',data%ocean_data%nzocn,zocn_dimid)
    else
       status = parallel_inq_dimid(NCO%id,'zocn',zocn_dimid)
    endif
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! Expanding restart variables: if 'restart' or 'hot' is present, we remove that
    ! word from the variable list, and flip the restartfile flag.  
    ! In CISM 2.0, 'restart' is the preferred name to represent restart variables, 
    ! but 'hot' is supported for backward compatibility.  Thus, we check for both.
    NCO%vars = ' '//trim(adjustl(NCO%vars))//' '  ! Need to maintain a space at beginning and end of list
    ! expanding restart variables
    pos = index(NCO%vars,' restart ') 
    if (pos.ne.0) then
       NCO%vars = NCO%vars(:pos)//NCO%vars(pos+8:)
       NCO%restartfile = .true.
    end if
    pos = index(NCO%vars,' hot ') 
    if (pos.ne.0) then
       NCO%vars = NCO%vars(:pos)//NCO%vars(pos+4:)
       NCO%restartfile = .true.
    end if
    ! Now apply necessary changes if the file is a restart file.
    if (NCO%restartfile) then
       if ((len_trim(NCO%vars) + len_trim(restart_variable_list) + 2) >= len(NCO%vars) ) then
          call write_log('Adding restart variables has made the list of output variables too long for file ' // NCO%filename, &
               GM_FATAL)
       else
          ! Expand the restart variable list 
          ! Need to maintain a space at beginning and end of list
          NCO%vars = trim(NCO%vars) // ' ' // trim(restart_variable_list) // ' ' ! (a module variable)  
          ! Set the xtype to be double (required for an exact restart)
          outfile%default_xtype = NF90_DOUBLE   
       endif
    end if

    ! Convert temp and flwa to versions on stag grid, if needed
    ! Note: this check must occur after restart variables are expanded which happens in glimmer_nc_readparams
    call check_for_tempstag(model%options%whichdycore,NCO)

    ! checking if we need to handle time averages
    pos = index(NCO%vars,"_tavg")
    if (pos.ne.0) then
       outfile%do_averages = .True.
    end if    

    ! Now that the output variable list is finalized, make sure we aren't truncating what the user intends to be output.
    ! Note: this only checks that the text in the variable list does not extend to within one character of the end of the variable.
    ! It does not handle the case where the user exactly fills the allowable length with variables or has a too-long list with more than one space between variable names.
    if ((len_trim(NCO%vars) + 1 ) >= len(NCO%vars)) then 
       call write_log('The list of output variables is too long for file ' // NCO%filename, GM_FATAL)
    endif


    ! MJH, 10/21/13: In the auto-generated code below, the creation of each output variable is wrapped by a check if the data for that 
    !   variable has a size greater than 0.  This is because of recently added checks in glide_types.F90 that don't fully allocate
    !   some variables if certain model options are disabled.  This is to lower memory requirements while running the model.
    !   The reason they have to be allocated with size zero rather than left unallocated is because the data for
    !   some netCDF output variables is defined with math, which causes an error if the operands are unallocated.
    !   Note that if a variable is not created, then it will not be subsequently written to.
    !   Also note that this change requires that data be a mandatory argument to this subroutine.

    ! Some output variables will need tavgf.  The value does not matter, but it must exist.  
    ! Nonetheless, for completeness give it the proper value that it has in glide_io_write.
    tavgf = outfile%total_time
    if (tavgf.ne.0.d0) then
       tavgf = 1.d0/tavgf
    end if
    ! Similarly, some output variables use the variable up.  Give it value of 0 here.
    up = 0

    !     level -- sigma layers
    if (.not.outfile%append) then
       call write_log('Creating variable level')
       status = parallel_def_var(NCO%id,'level',get_xtype(outfile,NF90_FLOAT), &
            (/level_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'sigma layers')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_sigma_coordinate')
       status = parallel_put_att(NCO%id, varid, 'positive', &
            'down')
     end if

    !     lithoz -- vertical coordinate of lithosphere layer
    if (.not.outfile%append) then
       call write_log('Creating variable lithoz')
       status = parallel_def_var(NCO%id,'lithoz',get_xtype(outfile,NF90_FLOAT), &
            (/lithoz_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertical coordinate of lithosphere layer')
     end if

    !     nlev_smb -- vertical level for SMB forcing
    if (.not.outfile%append) then
       call write_log('Creating variable nlev_smb')
       status = parallel_def_var(NCO%id,'nlev_smb',get_xtype(outfile,NF90_FLOAT), &
            (/nlev_smb_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertical level for SMB forcing')
     end if

    !     staglevel -- stag sigma layers
    if (.not.outfile%append) then
       call write_log('Creating variable staglevel')
       status = parallel_def_var(NCO%id,'staglevel',get_xtype(outfile,NF90_FLOAT), &
            (/staglevel_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'stag sigma layers')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_stag_sigma_coordinate')
       status = parallel_put_att(NCO%id, varid, 'positive', &
            'down')
     end if

    !     stagwbndlevel -- stag sigma layers with boundaries
    if (.not.outfile%append) then
       call write_log('Creating variable stagwbndlevel')
       status = parallel_def_var(NCO%id,'stagwbndlevel',get_xtype(outfile,NF90_FLOAT), &
            (/stagwbndlevel_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'stag sigma layers with boundaries')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_stag_sigma_coordinate_with_bnd')
       status = parallel_put_att(NCO%id, varid, 'positive', &
            'down')
     end if

    !     x0 -- Cartesian x-coordinate, velocity grid
    if (.not.outfile%append) then
       call write_log('Creating variable x0')
       status = parallel_def_var(NCO%id,'x0',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Cartesian x-coordinate, velocity grid')
       status = parallel_put_att(NCO%id, varid, 'axis', &
            'X')
     end if

    !     x1 -- Cartesian x-coordinate
    if (.not.outfile%append) then
       call write_log('Creating variable x1')
       status = parallel_def_var(NCO%id,'x1',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Cartesian x-coordinate')
       status = parallel_put_att(NCO%id, varid, 'axis', &
            'X')
     end if

    !     y0 -- Cartesian y-coordinate, velocity grid
    if (.not.outfile%append) then
       call write_log('Creating variable y0')
       status = parallel_def_var(NCO%id,'y0',get_xtype(outfile,NF90_FLOAT), &
            (/y0_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Cartesian y-coordinate, velocity grid')
       status = parallel_put_att(NCO%id, varid, 'axis', &
            'Y')
     end if

    !     y1 -- Cartesian y-coordinate
    if (.not.outfile%append) then
       call write_log('Creating variable y1')
       status = parallel_def_var(NCO%id,'y1',get_xtype(outfile,NF90_FLOAT), &
            (/y1_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Cartesian y-coordinate')
       status = parallel_put_att(NCO%id, varid, 'axis', &
            'Y')
     end if

    !     zocn -- ocean_z_coordinate
    if (.not.outfile%append) then
       call write_log('Creating variable zocn')
       status = parallel_def_var(NCO%id,'zocn',get_xtype(outfile,NF90_FLOAT), &
            (/zocn_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ocean_z_coordinate')
       status = parallel_put_att(NCO%id, varid, 'positive', &
            'up')
     end if

    !     S_ambient -- ambient ocean salinity
    pos = index(NCO%vars,' S_ambient ')
    status = parallel_inq_varid(NCO%id,'S_ambient',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%plume%S_ambient)) then
       call write_log('Creating variable S_ambient')
       status = parallel_def_var(NCO%id,'S_ambient',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'psu')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ambient ocean salinity')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'ambient ocean salinity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable S_ambient was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     T_ambient -- ambient ocean temperature
    pos = index(NCO%vars,' T_ambient ')
    status = parallel_inq_varid(NCO%id,'T_ambient',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%plume%T_ambient)) then
       call write_log('Creating variable T_ambient')
       status = parallel_def_var(NCO%id,'T_ambient',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ambient ocean temperature')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'ambient_ocean_temperature')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable T_ambient was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     acab -- accumulation, ablation rate
    pos = index(NCO%vars,' acab ')
    status = parallel_inq_varid(NCO%id,'acab',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%acab)) then
       call write_log('Creating variable acab')
       status = parallel_def_var(NCO%id,'acab',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year ice')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'accumulation, ablation rate')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable acab was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     acab_3d -- 3D surface mass balance
    pos = index(NCO%vars,' acab_3d ')
    status = parallel_inq_varid(NCO%id,'acab_3d',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%acab_3d)) then
       call write_log('Creating variable acab_3d')
       status = parallel_def_var(NCO%id,'acab_3d',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, nlev_smb_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'm/year ice')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            '3D surface mass balance')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_3d')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable acab_3d was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     acab_anomaly -- accumulation, ablation rate anomaly
    pos = index(NCO%vars,' acab_anomaly ')
    status = parallel_inq_varid(NCO%id,'acab_anomaly',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%acab_anomaly)) then
       call write_log('Creating variable acab_anomaly')
       status = parallel_def_var(NCO%id,'acab_anomaly',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'accumulation, ablation rate anomaly')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_anomaly')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable acab_anomaly was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     acab_applied -- applied accumulation, ablation rate
    pos = index(NCO%vars,' acab_applied ')
    status = parallel_inq_varid(NCO%id,'acab_applied',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%acab_applied)) then
       call write_log('Creating variable acab_applied')
       status = parallel_def_var(NCO%id,'acab_applied',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'applied accumulation, ablation rate')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_applied')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable acab_applied was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     acab_applied_tavg -- applied accumulation, ablation rate (time average)
    pos = index(NCO%vars,' acab_applied_tavg ')
    status = parallel_inq_varid(NCO%id,'acab_applied_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+17) = '                 '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%acab_applied_tavg)) then
       call write_log('Creating variable acab_applied_tavg')
       status = parallel_def_var(NCO%id,'acab_applied_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'applied accumulation, ablation rate (time average)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_applied')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean over years')
       status = parallel_put_att(NCO%id, varid, 'avg_factor', &
            'tavgf')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable acab_applied_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     acab_corrected -- corrected accumulation, ablation rate
    pos = index(NCO%vars,' acab_corrected ')
    status = parallel_inq_varid(NCO%id,'acab_corrected',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%acab_corrected)) then
       call write_log('Creating variable acab_corrected')
       status = parallel_def_var(NCO%id,'acab_corrected',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'corrected accumulation, ablation rate')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_corrected')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable acab_corrected was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     acab_gradz -- surface mass balance vertical gradient
    pos = index(NCO%vars,' acab_gradz ')
    status = parallel_inq_varid(NCO%id,'acab_gradz',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%acab_gradz)) then
       call write_log('Creating variable acab_gradz')
       status = parallel_def_var(NCO%id,'acab_gradz',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab/thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'm/year ice per m')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface mass balance vertical gradient')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_vertical_gradient')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable acab_gradz was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     acab_ref -- accumulation, ablation rate
    pos = index(NCO%vars,' acab_ref ')
    status = parallel_inq_varid(NCO%id,'acab_ref',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%acab_ref)) then
       call write_log('Creating variable acab_ref')
       status = parallel_def_var(NCO%id,'acab_ref',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year ice')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'accumulation, ablation rate')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_reference')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable acab_ref was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     adv_cfl_dt -- advective CFL maximum time step
    pos = index(NCO%vars,' adv_cfl_dt ')
    status = parallel_inq_varid(NCO%id,'adv_cfl_dt',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%numerics%adv_cfl_dt)) then
       call write_log('Creating variable adv_cfl_dt')
       status = parallel_def_var(NCO%id,'adv_cfl_dt',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'years')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'advective CFL maximum time step')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable adv_cfl_dt was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     area_factor -- area scale factor for stereographic map projection
    pos = index(NCO%vars,' area_factor ')
    status = parallel_inq_varid(NCO%id,'area_factor',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%projection%stere%area_factor)) then
       call write_log('Creating variable area_factor')
       status = parallel_def_var(NCO%id,'area_factor',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'area scale factor for stereographic map projection')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable area_factor was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     artm -- annual mean air temperature
    pos = index(NCO%vars,' artm ')
    status = parallel_inq_varid(NCO%id,'artm',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%artm)) then
       call write_log('Creating variable artm')
       status = parallel_def_var(NCO%id,'artm',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'annual mean air temperature')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'surface_air_temperature')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable artm was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     artm_3d -- 3D surface temperature
    pos = index(NCO%vars,' artm_3d ')
    status = parallel_inq_varid(NCO%id,'artm_3d',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%artm_3d)) then
       call write_log('Creating variable artm_3d')
       status = parallel_def_var(NCO%id,'artm_3d',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, nlev_smb_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'deg Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            '3D surface temperature')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_temperature_3d')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable artm_3d was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     artm_anomaly -- surface temperature anomaly
    pos = index(NCO%vars,' artm_anomaly ')
    status = parallel_inq_varid(NCO%id,'artm_anomaly',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%artm_anomaly)) then
       call write_log('Creating variable artm_anomaly')
       status = parallel_def_var(NCO%id,'artm_anomaly',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'deg Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface temperature anomaly')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_temperature_anomaly')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable artm_anomaly was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     artm_gradz -- surface temperature vertical gradient
    pos = index(NCO%vars,' artm_gradz ')
    status = parallel_inq_varid(NCO%id,'artm_gradz',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%artm_gradz)) then
       call write_log('Creating variable artm_gradz')
       status = parallel_def_var(NCO%id,'artm_gradz',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1./thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'deg Celsius per m')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface temperature vertical gradient')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_temperature_vertical_gradient')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable artm_gradz was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     artm_ref -- surface temperature at reference elevation
    pos = index(NCO%vars,' artm_ref ')
    status = parallel_inq_varid(NCO%id,'artm_ref',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%artm_ref)) then
       call write_log('Creating variable artm_ref')
       status = parallel_def_var(NCO%id,'artm_ref',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'deg Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface temperature at reference elevation')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_temperature_reference')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable artm_ref was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     basal_mbal_flux -- basal mass balance flux
    pos = index(NCO%vars,' basal_mbal_flux ')
    status = parallel_inq_varid(NCO%id,'basal_mbal_flux',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+15) = '               '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%basal_mbal_flux)) then
       call write_log('Creating variable basal_mbal_flux')
       status = parallel_def_var(NCO%id,'basal_mbal_flux',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/m2/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal mass balance flux')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_specific_mass_balance_flux')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable basal_mbal_flux was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     basal_mbal_flux_tavg -- basal mass balance flux (time average)
    pos = index(NCO%vars,' basal_mbal_flux_tavg ')
    status = parallel_inq_varid(NCO%id,'basal_mbal_flux_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+20) = '                    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%basal_mbal_flux_tavg)) then
       call write_log('Creating variable basal_mbal_flux_tavg')
       status = parallel_def_var(NCO%id,'basal_mbal_flux_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/m2/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal mass balance flux (time average)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_specific_mass_balance_flux')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean over years')
       status = parallel_put_att(NCO%id, varid, 'avg_factor', &
            'tavgf')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable basal_mbal_flux_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     basin_multiplier_array -- multiplier on a basin scale for the basal melting
    pos = index(NCO%vars,' basin_multiplier_array ')
    status = parallel_inq_varid(NCO%id,'basin_multiplier_array',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+22) = '                      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_melt%basin_multiplier_array)) then
       call write_log('Creating variable basin_multiplier_array')
       status = parallel_def_var(NCO%id,'basin_multiplier_array',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'multiplier on a basin scale for the basal melting')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable basin_multiplier_array was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     basin_number -- basin_number
    pos = index(NCO%vars,' basin_number ')
    status = parallel_inq_varid(NCO%id,'basin_number',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%ocean_data%basin_number)) then
       call write_log('Creating variable basin_number')
       status = parallel_def_var(NCO%id,'basin_number',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basin_number')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable basin_number was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     beta -- higher-order bed stress coefficient
    pos = index(NCO%vars,' beta ')
    status = parallel_inq_varid(NCO%id,'beta',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%beta)) then
       call write_log('Creating variable beta')
       status = parallel_def_var(NCO%id,'beta',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_beta))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa yr/m')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'higher-order bed stress coefficient')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable beta was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     beta_internal -- weighted higher-order bed stress coefficient
    pos = index(NCO%vars,' beta_internal ')
    status = parallel_inq_varid(NCO%id,'beta_internal',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+13) = '             '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%beta_internal)) then
       call write_log('Creating variable beta_internal')
       status = parallel_def_var(NCO%id,'beta_internal',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_beta))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa yr/m')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'weighted higher-order bed stress coefficient')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable beta_internal was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bfricflx -- basal friction heat flux
    pos = index(NCO%vars,' bfricflx ')
    status = parallel_inq_varid(NCO%id,'bfricflx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%bfricflx)) then
       call write_log('Creating variable bfricflx')
       status = parallel_def_var(NCO%id,'bfricflx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'watt/meter2')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal friction heat flux')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bfricflx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bheatflx -- downward basal heat flux
    pos = index(NCO%vars,' bheatflx ')
    status = parallel_inq_varid(NCO%id,'bheatflx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%bheatflx)) then
       call write_log('Creating variable bheatflx')
       status = parallel_def_var(NCO%id,'bheatflx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'watt/meter2')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'downward basal heat flux')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bheatflx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bmlt -- basal melt rate
    pos = index(NCO%vars,' bmlt ')
    status = parallel_inq_varid(NCO%id,'bmlt',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_melt%bmlt)) then
       call write_log('Creating variable bmlt')
       status = parallel_def_var(NCO%id,'bmlt',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal melt rate')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_melt_rate')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bmlt was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bmlt_applied -- applied basal melt rate
    pos = index(NCO%vars,' bmlt_applied ')
    status = parallel_inq_varid(NCO%id,'bmlt_applied',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_melt%bmlt_applied)) then
       call write_log('Creating variable bmlt_applied')
       status = parallel_def_var(NCO%id,'bmlt_applied',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'applied basal melt rate')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_melt_rate_applied')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bmlt_applied was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bmlt_applied_tavg -- applied basal melt rate (time average)
    pos = index(NCO%vars,' bmlt_applied_tavg ')
    status = parallel_inq_varid(NCO%id,'bmlt_applied_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+17) = '                 '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_melt%bmlt_applied_tavg)) then
       call write_log('Creating variable bmlt_applied_tavg')
       status = parallel_def_var(NCO%id,'bmlt_applied_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'applied basal melt rate (time average)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_melt_rate_applied')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean over years')
       status = parallel_put_att(NCO%id, varid, 'avg_factor', &
            'tavgf')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bmlt_applied_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bmlt_float -- basal melt rate for floating ice
    pos = index(NCO%vars,' bmlt_float ')
    status = parallel_inq_varid(NCO%id,'bmlt_float',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_melt%bmlt_float)) then
       call write_log('Creating variable bmlt_float')
       status = parallel_def_var(NCO%id,'bmlt_float',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal melt rate for floating ice')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'floating_ice_basal_melt_rate')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bmlt_float was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bmlt_float_anomaly -- basal melt rate anomaly for floating ice
    pos = index(NCO%vars,' bmlt_float_anomaly ')
    status = parallel_inq_varid(NCO%id,'bmlt_float_anomaly',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+18) = '                  '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_melt%bmlt_float_anomaly)) then
       call write_log('Creating variable bmlt_float_anomaly')
       status = parallel_def_var(NCO%id,'bmlt_float_anomaly',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal melt rate anomaly for floating ice')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'floating_ice_basal_melt_rate_anomaly')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bmlt_float_anomaly was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bmlt_float_external -- external basal melt rate for floating ice
    pos = index(NCO%vars,' bmlt_float_external ')
    status = parallel_inq_varid(NCO%id,'bmlt_float_external',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+19) = '                   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_melt%bmlt_float_external)) then
       call write_log('Creating variable bmlt_float_external')
       status = parallel_def_var(NCO%id,'bmlt_float_external',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'external basal melt rate for floating ice')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'floating_ice_basal_melt_rate_external')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bmlt_float_external was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bmlt_ground -- basal melt rate for grounded ice
    pos = index(NCO%vars,' bmlt_ground ')
    status = parallel_inq_varid(NCO%id,'bmlt_ground',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_melt%bmlt_ground)) then
       call write_log('Creating variable bmlt_ground')
       status = parallel_def_var(NCO%id,'bmlt_ground',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal melt rate for grounded ice')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'grounded_ice_basal_melt_rate')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bmlt_ground was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bpmp -- basal pressure melting point temperature
    pos = index(NCO%vars,' bpmp ')
    status = parallel_inq_varid(NCO%id,'bpmp',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%bpmp)) then
       call write_log('Creating variable bpmp')
       status = parallel_def_var(NCO%id,'bpmp',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal pressure melting point temperature')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'pressure_melting_point_temperature')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bpmp was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btemp -- basal ice temperature
    pos = index(NCO%vars,' btemp ')
    status = parallel_inq_varid(NCO%id,'btemp',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%btemp)) then
       call write_log('Creating variable btemp')
       status = parallel_def_var(NCO%id,'btemp',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal ice temperature')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_temperature')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btemp was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btemp_float -- basal ice temperature for floating ice
    pos = index(NCO%vars,' btemp_float ')
    status = parallel_inq_varid(NCO%id,'btemp_float',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%btemp_float)) then
       call write_log('Creating variable btemp_float')
       status = parallel_def_var(NCO%id,'btemp_float',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal ice temperature for floating ice')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_temperature_float')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btemp_float was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btemp_ground -- basal ice temperature for grounded ice
    pos = index(NCO%vars,' btemp_ground ')
    status = parallel_inq_varid(NCO%id,'btemp_ground',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%btemp_ground)) then
       call write_log('Creating variable btemp_ground')
       status = parallel_def_var(NCO%id,'btemp_ground',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal ice temperature for grounded ice')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_temperature_ground')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btemp_ground was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btract -- basal traction (magnitude)
    pos = index(NCO%vars,' btract ')
    status = parallel_inq_varid(NCO%id,'btract',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%btract)) then
       call write_log('Creating variable btract')
       status = parallel_def_var(NCO%id,'btract',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal traction (magnitude)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btract was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btractx -- basal traction (x-direction comp)
    pos = index(NCO%vars,' btractx ')
    status = parallel_inq_varid(NCO%id,'btractx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%btractx)) then
       call write_log('Creating variable btractx')
       status = parallel_def_var(NCO%id,'btractx',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal traction (x-direction comp)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btractx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btractx_extend -- basal traction (x-direction comp)
    pos = index(NCO%vars,' btractx_extend ')
    status = parallel_inq_varid(NCO%id,'btractx_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%btractx_extend)) then
       call write_log('Creating variable btractx_extend')
       status = parallel_def_var(NCO%id,'btractx_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal traction (x-direction comp)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btractx_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btracty -- basal traction (y-direction comp)
    pos = index(NCO%vars,' btracty ')
    status = parallel_inq_varid(NCO%id,'btracty',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%btracty)) then
       call write_log('Creating variable btracty')
       status = parallel_def_var(NCO%id,'btracty',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal traction (y-direction comp)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btracty was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btracty_extend -- basal traction (y-direction comp)
    pos = index(NCO%vars,' btracty_extend ')
    status = parallel_inq_varid(NCO%id,'btracty_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%btracty_extend)) then
       call write_log('Creating variable btracty_extend')
       status = parallel_def_var(NCO%id,'btracty_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal traction (y-direction comp)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btracty_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     btrc -- basal slip coefficient
    pos = index(NCO%vars,' btrc ')
    status = parallel_inq_varid(NCO%id,'btrc',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%btrc)) then
       call write_log('Creating variable btrc')
       status = parallel_def_var(NCO%id,'btrc',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_btrc))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/pascal/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal slip coefficient')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable btrc was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bwat -- basal water depth
    pos = index(NCO%vars,' bwat ')
    status = parallel_inq_varid(NCO%id,'bwat',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_hydro%bwat)) then
       call write_log('Creating variable bwat')
       status = parallel_def_var(NCO%id,'bwat',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal water depth')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bwat was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bwat_diag -- basal water depth, that does not influence the temperature, so it can work with bwatflx routing
    pos = index(NCO%vars,' bwat_diag ')
    status = parallel_inq_varid(NCO%id,'bwat_diag',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_hydro%bwat_diag)) then
       call write_log('Creating variable bwat_diag')
       status = parallel_def_var(NCO%id,'bwat_diag',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal water depth, that does not influence the temperature, so it can work with bwatflx routing')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bwat_diag was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     bwatflx -- basal water flux
    pos = index(NCO%vars,' bwatflx ')
    status = parallel_inq_varid(NCO%id,'bwatflx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_hydro%bwatflx)) then
       call write_log('Creating variable bwatflx')
       status = parallel_def_var(NCO%id,'bwatflx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal water flux')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable bwatflx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     c_flux_array -- spatially variant flux to depth constant to allow turbulent flow
    pos = index(NCO%vars,' c_flux_array ')
    status = parallel_inq_varid(NCO%id,'c_flux_array',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_hydro%c_flux_array)) then
       call write_log('Creating variable c_flux_array')
       status = parallel_def_var(NCO%id,'c_flux_array',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '-')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'spatially variant flux to depth constant to allow turbulent flow')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable c_flux_array was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     c_space_factor -- spatial factor for basal shear stress
    pos = index(NCO%vars,' c_space_factor ')
    status = parallel_inq_varid(NCO%id,'c_space_factor',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%c_space_factor)) then
       call write_log('Creating variable c_space_factor')
       status = parallel_def_var(NCO%id,'c_space_factor',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'spatial factor for basal shear stress')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable c_space_factor was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     calving_flux -- calving flux
    pos = index(NCO%vars,' calving_flux ')
    status = parallel_inq_varid(NCO%id,'calving_flux',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%calving_flux)) then
       call write_log('Creating variable calving_flux')
       status = parallel_def_var(NCO%id,'calving_flux',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/m2/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'calving flux')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_specific_mass_flux_due_to_calving')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable calving_flux was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     calving_flux_tavg -- calving flux (time average)
    pos = index(NCO%vars,' calving_flux_tavg ')
    status = parallel_inq_varid(NCO%id,'calving_flux_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+17) = '                 '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%calving_flux_tavg)) then
       call write_log('Creating variable calving_flux_tavg')
       status = parallel_def_var(NCO%id,'calving_flux_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/m2/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'calving flux (time average)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_specific_mass_flux_due_to_calving')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean over years')
       status = parallel_put_att(NCO%id, varid, 'avg_factor', &
            'tavgf')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable calving_flux_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     calving_lateral -- lateral calving rate
    pos = index(NCO%vars,' calving_lateral ')
    status = parallel_inq_varid(NCO%id,'calving_lateral',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+15) = '               '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%lateral_rate)) then
       call write_log('Creating variable calving_lateral')
       status = parallel_def_var(NCO%id,'calving_lateral',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'lateral calving rate')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable calving_lateral was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     calving_mask -- calving mask
    pos = index(NCO%vars,' calving_mask ')
    status = parallel_inq_varid(NCO%id,'calving_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%calving_mask)) then
       call write_log('Creating variable calving_mask')
       status = parallel_def_var(NCO%id,'calving_mask',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'calving mask')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable calving_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     calving_rate -- rate of mass loss by calving
    pos = index(NCO%vars,' calving_rate ')
    status = parallel_inq_varid(NCO%id,'calving_rate',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%calving_rate)) then
       call write_log('Creating variable calving_rate')
       status = parallel_def_var(NCO%id,'calving_rate',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'rate of mass loss by calving')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable calving_rate was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     calving_rate_tavg -- rate of mass loss by calving (time average)
    pos = index(NCO%vars,' calving_rate_tavg ')
    status = parallel_inq_varid(NCO%id,'calving_rate_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+17) = '                 '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%calving_rate_tavg)) then
       call write_log('Creating variable calving_rate_tavg')
       status = parallel_def_var(NCO%id,'calving_rate_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'rate of mass loss by calving (time average)')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean over years')
       status = parallel_put_att(NCO%id, varid, 'avg_factor', &
            'tavgf')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable calving_rate_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     calving_thck -- thickness of calving ice
    pos = index(NCO%vars,' calving_thck ')
    status = parallel_inq_varid(NCO%id,'calving_thck',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%calving_thck)) then
       call write_log('Creating variable calving_thck')
       status = parallel_def_var(NCO%id,'calving_thck',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'thickness of calving ice')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable calving_thck was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     cell_area -- cell area of cism grid
    pos = index(NCO%vars,' cell_area ')
    status = parallel_inq_varid(NCO%id,'cell_area',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%cell_area)) then
       call write_log('Creating variable cell_area')
       status = parallel_def_var(NCO%id,'cell_area',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(len0*len0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter2')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'cell area of cism grid')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable cell_area was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     coulomb_c -- spatially varying C for Coulomb sliding, staggered grid
    pos = index(NCO%vars,' coulomb_c ')
    status = parallel_inq_varid(NCO%id,'coulomb_c',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%coulomb_c)) then
       call write_log('Creating variable coulomb_c')
       status = parallel_def_var(NCO%id,'coulomb_c',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'spatially varying C for Coulomb sliding, staggered grid')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable coulomb_c was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     coulomb_c_relax -- target for coulomb c to relax towards
    pos = index(NCO%vars,' coulomb_c_relax ')
    status = parallel_inq_varid(NCO%id,'coulomb_c_relax',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+15) = '               '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%coulomb_c_relax)) then
       call write_log('Creating variable coulomb_c_relax')
       status = parallel_def_var(NCO%id,'coulomb_c_relax',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'target for coulomb c to relax towards')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable coulomb_c_relax was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     damage -- ice damage
    pos = index(NCO%vars,' damage ')
    status = parallel_inq_varid(NCO%id,'damage',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%damage)) then
       call write_log('Creating variable damage')
       status = parallel_def_var(NCO%id,'damage',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless [0,1]')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice damage')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable damage was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     deltaT_ocn -- deltaT_ocn
    pos = index(NCO%vars,' deltaT_ocn ')
    status = parallel_inq_varid(NCO%id,'deltaT_ocn',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%ocean_data%deltaT_ocn)) then
       call write_log('Creating variable deltaT_ocn')
       status = parallel_def_var(NCO%id,'deltaT_ocn',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degrees K')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'deltaT_ocn')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable deltaT_ocn was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     diff_cfl_dt -- diffusive CFL maximum time step
    pos = index(NCO%vars,' diff_cfl_dt ')
    status = parallel_inq_varid(NCO%id,'diff_cfl_dt',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%numerics%diff_cfl_dt)) then
       call write_log('Creating variable diff_cfl_dt')
       status = parallel_def_var(NCO%id,'diff_cfl_dt',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'years')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'diffusive CFL maximum time step')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable diff_cfl_dt was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     diffu -- apparent diffusivity
    pos = index(NCO%vars,' diffu ')
    status = parallel_inq_varid(NCO%id,'diffu',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%diffu)) then
       call write_log('Creating variable diffu')
       status = parallel_def_var(NCO%id,'diffu',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_diffu))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter2/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'apparent diffusivity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable diffu was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dissip -- dissipation rate (W m-3) divided by rhoi Ci
    pos = index(NCO%vars,' dissip ')
    status = parallel_inq_varid(NCO%id,'dissip',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%dissip)) then
       call write_log('Creating variable dissip')
       status = parallel_def_var(NCO%id,'dissip',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'deg C/yr')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'dissipation rate (W m-3) divided by rhoi Ci')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable dissip was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dissipstag -- dissipation rate (W m-3) divided by rhoi Ci
    pos = index(NCO%vars,' dissipstag ')
    status = parallel_inq_varid(NCO%id,'dissipstag',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%dissip)) then
       call write_log('Creating variable dissipstag')
       status = parallel_def_var(NCO%id,'dissipstag',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'deg C/yr')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'dissipation rate (W m-3) divided by rhoi Ci')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable dissipstag was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     divu -- horizontal divergence
    pos = index(NCO%vars,' divu ')
    status = parallel_inq_varid(NCO%id,'divu',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%divu)) then
       call write_log('Creating variable divu')
       status = parallel_def_var(NCO%id,'divu',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'horizontal divergence')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable divu was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dthck_dt -- rate of ice thickness change
    pos = index(NCO%vars,' dthck_dt ')
    status = parallel_inq_varid(NCO%id,'dthck_dt',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%dthck_dt)) then
       call write_log('Creating variable dthck_dt')
       status = parallel_def_var(NCO%id,'dthck_dt',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'rate of ice thickness change')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable dthck_dt was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dthck_dt_obs -- observed rate of ice thickness change
    pos = index(NCO%vars,' dthck_dt_obs ')
    status = parallel_inq_varid(NCO%id,'dthck_dt_obs',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%dthck_dt_obs)) then
       call write_log('Creating variable dthck_dt_obs')
       status = parallel_def_var(NCO%id,'dthck_dt_obs',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'observed rate of ice thickness change')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable dthck_dt_obs was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dthck_dt_obs_basin -- observed rate of ice thickness change, basin average
    pos = index(NCO%vars,' dthck_dt_obs_basin ')
    status = parallel_inq_varid(NCO%id,'dthck_dt_obs_basin',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+18) = '                  '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%dthck_dt_obs_basin)) then
       call write_log('Creating variable dthck_dt_obs_basin')
       status = parallel_def_var(NCO%id,'dthck_dt_obs_basin',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'observed rate of ice thickness change, basin average')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable dthck_dt_obs_basin was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dthckdtm -- tendency of ice thickness (NOTE: Glide only)
    pos = index(NCO%vars,' dthckdtm ')
    status = parallel_inq_varid(NCO%id,'dthckdtm',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geomderv%dthckdtm)) then
       call write_log('Creating variable dthckdtm')
       status = parallel_def_var(NCO%id,'dthckdtm',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'tendency of ice thickness (NOTE: Glide only)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable dthckdtm was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     dusrfdtm -- rate of upper ice surface elevation change (NOTE: Glide only)
    pos = index(NCO%vars,' dusrfdtm ')
    status = parallel_inq_varid(NCO%id,'dusrfdtm',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geomderv%dusrfdtm)) then
       call write_log('Creating variable dusrfdtm')
       status = parallel_def_var(NCO%id,'dusrfdtm',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_acab))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'rate of upper ice surface elevation change (NOTE: Glide only)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable dusrfdtm was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     effecpress -- effective pressure
    pos = index(NCO%vars,' effecpress ')
    status = parallel_inq_varid(NCO%id,'effecpress',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%effecpress)) then
       call write_log('Creating variable effecpress')
       status = parallel_def_var(NCO%id,'effecpress',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'effective pressure')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable effecpress was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     efvs -- effective viscosity
    pos = index(NCO%vars,' efvs ')
    status = parallel_inq_varid(NCO%id,'efvs',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%efvs)) then
       call write_log('Creating variable efvs')
       status = parallel_def_var(NCO%id,'efvs',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_efvs))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pascal * years')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'effective viscosity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable efvs was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     enthalpy -- specific enthalpy
    pos = index(NCO%vars,' enthalpy ')
    status = parallel_inq_varid(NCO%id,'enthalpy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%enthalpy)) then
       call write_log('Creating variable enthalpy')
       status = parallel_def_var(NCO%id,'enthalpy',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, stagwbndlevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'J/m3')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'specific enthalpy')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable enthalpy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     eps_eff -- effective strain rate
    pos = index(NCO%vars,' eps_eff ')
    status = parallel_inq_varid(NCO%id,'eps_eff',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%strain_rate%scalar)) then
       call write_log('Creating variable eps_eff')
       status = parallel_def_var(NCO%id,'eps_eff',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'effective strain rate')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable eps_eff was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     eps_eigen1 -- first eigenvalue of horizontal strain rate tensor
    pos = index(NCO%vars,' eps_eigen1 ')
    status = parallel_inq_varid(NCO%id,'eps_eigen1',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%eps_eigen1)) then
       call write_log('Creating variable eps_eigen1')
       status = parallel_def_var(NCO%id,'eps_eigen1',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'first eigenvalue of horizontal strain rate tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable eps_eigen1 was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     eps_eigen2 -- second eigenvalue of horizontal strain rate tensor
    pos = index(NCO%vars,' eps_eigen2 ')
    status = parallel_inq_varid(NCO%id,'eps_eigen2',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%eps_eigen2)) then
       call write_log('Creating variable eps_eigen2')
       status = parallel_def_var(NCO%id,'eps_eigen2',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'second eigenvalue of horizontal strain rate tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable eps_eigen2 was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     eps_xx -- xx component of strain rate tensor
    pos = index(NCO%vars,' eps_xx ')
    status = parallel_inq_varid(NCO%id,'eps_xx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%strain_rate%xx)) then
       call write_log('Creating variable eps_xx')
       status = parallel_def_var(NCO%id,'eps_xx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'xx component of strain rate tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable eps_xx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     eps_xy -- xy component of strain rate tensor
    pos = index(NCO%vars,' eps_xy ')
    status = parallel_inq_varid(NCO%id,'eps_xy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%strain_rate%xy)) then
       call write_log('Creating variable eps_xy')
       status = parallel_def_var(NCO%id,'eps_xy',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'xy component of strain rate tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable eps_xy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     eps_xz -- xz component of strain rate tensor
    pos = index(NCO%vars,' eps_xz ')
    status = parallel_inq_varid(NCO%id,'eps_xz',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%strain_rate%xz)) then
       call write_log('Creating variable eps_xz')
       status = parallel_def_var(NCO%id,'eps_xz',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'xz component of strain rate tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable eps_xz was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     eps_yy -- yy component of strain rate tensor
    pos = index(NCO%vars,' eps_yy ')
    status = parallel_inq_varid(NCO%id,'eps_yy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%strain_rate%yy)) then
       call write_log('Creating variable eps_yy')
       status = parallel_def_var(NCO%id,'eps_yy',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'yy component of strain rate tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable eps_yy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     eps_yz -- yz component of strain rate tensor
    pos = index(NCO%vars,' eps_yz ')
    status = parallel_inq_varid(NCO%id,'eps_yz',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%strain_rate%yz)) then
       call write_log('Creating variable eps_yz')
       status = parallel_def_var(NCO%id,'eps_yz',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'yz component of strain rate tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable eps_yz was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     eus -- global average sea level
    pos = index(NCO%vars,' eus ')
    status = parallel_inq_varid(NCO%id,'eus',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+3) = '   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%eus)) then
       call write_log('Creating variable eus')
       status = parallel_def_var(NCO%id,'eus',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'global average sea level')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'global_average_sea_level_change')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable eus was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     f_effecpress_bwat -- effective pressure factor from bwatflx
    pos = index(NCO%vars,' f_effecpress_bwat ')
    status = parallel_inq_varid(NCO%id,'f_effecpress_bwat',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+17) = '                 '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%f_effecpress_bwat)) then
       call write_log('Creating variable f_effecpress_bwat')
       status = parallel_def_var(NCO%id,'f_effecpress_bwat',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'effective pressure factor from bwatflx')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable f_effecpress_bwat was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     f_effecpress_bwat_target -- Relaxation target of effective pressure ratio
    pos = index(NCO%vars,' f_effecpress_bwat_target ')
    status = parallel_inq_varid(NCO%id,'f_effecpress_bwat_target',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+24) = '                        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%f_effecpress_bwat_target)) then
       call write_log('Creating variable f_effecpress_bwat_target')
       status = parallel_def_var(NCO%id,'f_effecpress_bwat_target',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Relaxation target of effective pressure ratio')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable f_effecpress_bwat_target was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     f_effecpress_bwatflx -- effective pressure factor from the bwatflx
    pos = index(NCO%vars,' f_effecpress_bwatflx ')
    status = parallel_inq_varid(NCO%id,'f_effecpress_bwatflx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+20) = '                    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%f_effecpress_bwatflx)) then
       call write_log('Creating variable f_effecpress_bwatflx')
       status = parallel_def_var(NCO%id,'f_effecpress_bwatflx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'effective pressure factor from the bwatflx')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable f_effecpress_bwatflx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     f_effecpress_ocean_p -- effective pressure factor from ocean_p
    pos = index(NCO%vars,' f_effecpress_ocean_p ')
    status = parallel_inq_varid(NCO%id,'f_effecpress_ocean_p',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+20) = '                    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%f_effecpress_ocean_p)) then
       call write_log('Creating variable f_effecpress_ocean_p')
       status = parallel_def_var(NCO%id,'f_effecpress_ocean_p',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'effective pressure factor from ocean_p')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable f_effecpress_ocean_p was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     f_flotation -- flotation function
    pos = index(NCO%vars,' f_flotation ')
    status = parallel_inq_varid(NCO%id,'f_flotation',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%f_flotation)) then
       call write_log('Creating variable f_flotation')
       status = parallel_def_var(NCO%id,'f_flotation',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'flotation function')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable f_flotation was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     f_ground -- grounded ice fraction
    pos = index(NCO%vars,' f_ground ')
    status = parallel_inq_varid(NCO%id,'f_ground',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%f_ground)) then
       call write_log('Creating variable f_ground')
       status = parallel_def_var(NCO%id,'f_ground',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless [0,1]')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'grounded ice fraction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'grounded_fraction')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable f_ground was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     f_ground_cell -- grounded ice fraction in cells
    pos = index(NCO%vars,' f_ground_cell ')
    status = parallel_inq_varid(NCO%id,'f_ground_cell',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+13) = '             '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%f_ground_cell)) then
       call write_log('Creating variable f_ground_cell')
       status = parallel_def_var(NCO%id,'f_ground_cell',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless [0,1]')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'grounded ice fraction in cells')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'grounded_fraction_cell')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable f_ground_cell was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     f_ground_obs -- grounded ice fraction
    pos = index(NCO%vars,' f_ground_obs ')
    status = parallel_inq_varid(NCO%id,'f_ground_obs',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%f_ground_obs)) then
       call write_log('Creating variable f_ground_obs')
       status = parallel_def_var(NCO%id,'f_ground_obs',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless [0,1]')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'grounded ice fraction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'grounded_fraction')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable f_ground_obs was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ff_invert_mask -- mask where to increase the flow factor to simulate damage
    pos = index(NCO%vars,' ff_invert_mask ')
    status = parallel_inq_varid(NCO%id,'ff_invert_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%inversion%ff_invert_mask)) then
       call write_log('Creating variable ff_invert_mask')
       status = parallel_def_var(NCO%id,'ff_invert_mask',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '-')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask where to increase the flow factor to simulate damage')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ff_invert_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     floating_mask -- mask for floating ice
    pos = index(NCO%vars,' floating_mask ')
    status = parallel_inq_varid(NCO%id,'floating_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+13) = '             '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%floating_mask)) then
       call write_log('Creating variable floating_mask')
       status = parallel_def_var(NCO%id,'floating_mask',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask for floating ice')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable floating_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     floating_thck_target -- target thickness for floating ice
    pos = index(NCO%vars,' floating_thck_target ')
    status = parallel_inq_varid(NCO%id,'floating_thck_target',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+20) = '                    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%inversion%floating_thck_target)) then
       call write_log('Creating variable floating_thck_target')
       status = parallel_def_var(NCO%id,'floating_thck_target',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'target thickness for floating ice')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable floating_thck_target was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     flow_enhancement_factor -- flow enhancement factor
    pos = index(NCO%vars,' flow_enhancement_factor ')
    status = parallel_inq_varid(NCO%id,'flow_enhancement_factor',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+23) = '                       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%flow_enhancement_factor)) then
       call write_log('Creating variable flow_enhancement_factor')
       status = parallel_def_var(NCO%id,'flow_enhancement_factor',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'flow enhancement factor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable flow_enhancement_factor was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     flwa -- Pre-exponential flow law parameter
    pos = index(NCO%vars,' flwa ')
    status = parallel_inq_varid(NCO%id,'flwa',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%flwa)) then
       call write_log('Creating variable flwa')
       status = parallel_def_var(NCO%id,'flwa',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_flwa))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'pascal**(-n) year**(-1)')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Pre-exponential flow law parameter')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable flwa was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     flwastag -- Pre-exponential flow law parameter
    pos = index(NCO%vars,' flwastag ')
    status = parallel_inq_varid(NCO%id,'flwastag',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%flwa)) then
       call write_log('Creating variable flwastag')
       status = parallel_def_var(NCO%id,'flwastag',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_flwa))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'pascal**(-n) year**(-1)')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Pre-exponential flow law parameter')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable flwastag was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     gl_flux -- grounding line flux
    pos = index(NCO%vars,' gl_flux ')
    status = parallel_inq_varid(NCO%id,'gl_flux',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%gl_flux)) then
       call write_log('Creating variable gl_flux')
       status = parallel_def_var(NCO%id,'gl_flux',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/m/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'grounding line flux')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_mass_flux_at_grounding_line')
       status = parallel_put_att(NCO%id, varid, 'coordinate', &
            'lon lat')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable gl_flux was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     gl_flux_east -- grounding line flux eastward
    pos = index(NCO%vars,' gl_flux_east ')
    status = parallel_inq_varid(NCO%id,'gl_flux_east',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%gl_flux_east)) then
       call write_log('Creating variable gl_flux_east')
       status = parallel_def_var(NCO%id,'gl_flux_east',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/m/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'grounding line flux eastward')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_mass_flux_at_grounding_line_eastward')
       status = parallel_put_att(NCO%id, varid, 'coordinate', &
            'lon lat')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable gl_flux_east was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     gl_flux_north -- grounding line flux northward
    pos = index(NCO%vars,' gl_flux_north ')
    status = parallel_inq_varid(NCO%id,'gl_flux_north',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+13) = '             '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%gl_flux_north)) then
       call write_log('Creating variable gl_flux_north')
       status = parallel_def_var(NCO%id,'gl_flux_north',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/m/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'grounding line flux northward')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_mass_flux_at_grounding_line_northward')
       status = parallel_put_att(NCO%id, varid, 'coordinate', &
            'lon lat')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable gl_flux_north was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     gl_flux_tavg -- grounding line flux (time average)
    pos = index(NCO%vars,' gl_flux_tavg ')
    status = parallel_inq_varid(NCO%id,'gl_flux_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%gl_flux_tavg)) then
       call write_log('Creating variable gl_flux_tavg')
       status = parallel_def_var(NCO%id,'gl_flux_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/m/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'grounding line flux (time average)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_mass_flux_at_grounding_line')
       status = parallel_put_att(NCO%id, varid, 'coordinate', &
            'lon lat')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean over years')
       status = parallel_put_att(NCO%id, varid, 'avg_factor', &
            'tavgf')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable gl_flux_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     gravity -- gravitational acceleration
    pos = index(NCO%vars,' gravity ')
    status = parallel_inq_varid(NCO%id,'gravity',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(grav)) then
       call write_log('Creating variable gravity')
       status = parallel_def_var(NCO%id,'gravity',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/s/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'gravitational acceleration')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'gravity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable gravity was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     grounded_mask -- mask for grounded ice
    pos = index(NCO%vars,' grounded_mask ')
    status = parallel_inq_varid(NCO%id,'grounded_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+13) = '             '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%grounded_mask)) then
       call write_log('Creating variable grounded_mask')
       status = parallel_def_var(NCO%id,'grounded_mask',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask for grounded ice')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable grounded_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     head -- hydraulic head
    pos = index(NCO%vars,' head ')
    status = parallel_inq_varid(NCO%id,'head',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_hydro%head)) then
       call write_log('Creating variable head')
       status = parallel_def_var(NCO%id,'head',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'hydraulic head')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable head was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     iarea -- area covered by ice
    pos = index(NCO%vars,' iarea ')
    status = parallel_inq_varid(NCO%id,'iarea',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%iarea)) then
       call write_log('Creating variable iarea')
       status = parallel_def_var(NCO%id,'iarea',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'm2')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'area covered by ice')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable iarea was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     iareaf -- area covered by floating ice
    pos = index(NCO%vars,' iareaf ')
    status = parallel_inq_varid(NCO%id,'iareaf',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%iareaf)) then
       call write_log('Creating variable iareaf')
       status = parallel_def_var(NCO%id,'iareaf',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'm2')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'area covered by floating ice')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable iareaf was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     iareag -- area covered by grounded ice
    pos = index(NCO%vars,' iareag ')
    status = parallel_inq_varid(NCO%id,'iareag',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%iareag)) then
       call write_log('Creating variable iareag')
       status = parallel_def_var(NCO%id,'iareag',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'm2')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'area covered by grounded ice')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable iareag was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_age -- ice age
    pos = index(NCO%vars,' ice_age ')
    status = parallel_inq_varid(NCO%id,'ice_age',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%ice_age)) then
       call write_log('Creating variable ice_age')
       status = parallel_def_var(NCO%id,'ice_age',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(tim0/scyr))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice age')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_age')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ice_age was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_cap_mask -- mask for ice caps
    pos = index(NCO%vars,' ice_cap_mask ')
    status = parallel_inq_varid(NCO%id,'ice_cap_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%ice_cap_mask)) then
       call write_log('Creating variable ice_cap_mask')
       status = parallel_def_var(NCO%id,'ice_cap_mask',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask for ice caps')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ice_cap_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_domain_mask -- mask for where ice is potentially present and active
    pos = index(NCO%vars,' ice_domain_mask ')
    status = parallel_inq_varid(NCO%id,'ice_domain_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+15) = '               '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%general%ice_domain_mask)) then
       call write_log('Creating variable ice_domain_mask')
       status = parallel_def_var(NCO%id,'ice_domain_mask',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask for where ice is potentially present and active')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ice_domain_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_fraction_retreat_mask -- real mask for forced ice retreat
    pos = index(NCO%vars,' ice_fraction_retreat_mask ')
    status = parallel_inq_varid(NCO%id,'ice_fraction_retreat_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+25) = '                         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%ice_fraction_retreat_mask)) then
       call write_log('Creating variable ice_fraction_retreat_mask')
       status = parallel_def_var(NCO%id,'ice_fraction_retreat_mask',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'real mask for forced ice retreat')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ice_fraction_retreat_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_mask -- mask for ice (1) or no ice (0)
    pos = index(NCO%vars,' ice_mask ')
    status = parallel_inq_varid(NCO%id,'ice_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%ice_mask)) then
       call write_log('Creating variable ice_mask')
       status = parallel_def_var(NCO%id,'ice_mask',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask for ice (1) or no ice (0)')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ice_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_mask_stag -- mask on staggered grid for ice (1) or no ice (0)
    pos = index(NCO%vars,' ice_mask_stag ')
    status = parallel_inq_varid(NCO%id,'ice_mask_stag',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+13) = '             '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%ice_mask_stag)) then
       call write_log('Creating variable ice_mask_stag')
       status = parallel_def_var(NCO%id,'ice_mask_stag',get_xtype(outfile,NF90_INT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask on staggered grid for ice (1) or no ice (0)')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ice_mask_stag was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_sheet_mask -- mask for ice sheet
    pos = index(NCO%vars,' ice_sheet_mask ')
    status = parallel_inq_varid(NCO%id,'ice_sheet_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%ice_sheet_mask)) then
       call write_log('Creating variable ice_sheet_mask')
       status = parallel_def_var(NCO%id,'ice_sheet_mask',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask for ice sheet')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ice_sheet_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_specific_heat -- ice specific heat
    pos = index(NCO%vars,' ice_specific_heat ')
    status = parallel_inq_varid(NCO%id,'ice_specific_heat',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+17) = '                 '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(shci)) then
       call write_log('Creating variable ice_specific_heat')
       status = parallel_def_var(NCO%id,'ice_specific_heat',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'J/kg/K')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice specific heat')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'ice_specific_heat')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ice_specific_heat was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ice_thermal_conductivity -- ice thermal conductivity
    pos = index(NCO%vars,' ice_thermal_conductivity ')
    status = parallel_inq_varid(NCO%id,'ice_thermal_conductivity',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+24) = '                        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(coni)) then
       call write_log('Creating variable ice_thermal_conductivity')
       status = parallel_def_var(NCO%id,'ice_thermal_conductivity',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'J/(K kg)')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice thermal conductivity')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'ice_thermal_conductivity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ice_thermal_conductivity was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     imass -- ice mass
    pos = index(NCO%vars,' imass ')
    status = parallel_inq_varid(NCO%id,'imass',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%imass)) then
       call write_log('Creating variable imass')
       status = parallel_def_var(NCO%id,'imass',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice mass')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable imass was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     imass_above_flotation -- ice mass above flotation
    pos = index(NCO%vars,' imass_above_flotation ')
    status = parallel_inq_varid(NCO%id,'imass_above_flotation',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+21) = '                     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%imass_above_flotation)) then
       call write_log('Creating variable imass_above_flotation')
       status = parallel_def_var(NCO%id,'imass_above_flotation',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice mass above flotation')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable imass_above_flotation was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ivol -- ice volume
    pos = index(NCO%vars,' ivol ')
    status = parallel_inq_varid(NCO%id,'ivol',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%ivol)) then
       call write_log('Creating variable ivol')
       status = parallel_def_var(NCO%id,'ivol',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'm3')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice volume')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ivol was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     kinbcmask -- Mask of locations where uvel, vvel value should be held constant
    pos = index(NCO%vars,' kinbcmask ')
    status = parallel_inq_varid(NCO%id,'kinbcmask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%kinbcmask)) then
       call write_log('Creating variable kinbcmask')
       status = parallel_def_var(NCO%id,'kinbcmask',get_xtype(outfile,NF90_INT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Mask of locations where uvel, vvel value should be held constant')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable kinbcmask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     lat -- latitude
    pos = index(NCO%vars,' lat ')
    status = parallel_inq_varid(NCO%id,'lat',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+3) = '   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%general%lat)) then
       call write_log('Creating variable lat')
       status = parallel_def_var(NCO%id,'lat',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degreeN')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'latitude')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'latitude')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable lat was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     litho_temp -- lithosphere temperature
    pos = index(NCO%vars,' litho_temp ')
    status = parallel_inq_varid(NCO%id,'litho_temp',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%lithot%temp)) then
       call write_log('Creating variable litho_temp')
       status = parallel_def_var(NCO%id,'litho_temp',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, lithoz_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'lithosphere temperature')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable litho_temp was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     load -- bedrock deflection from applied load
    pos = index(NCO%vars,' load ')
    status = parallel_inq_varid(NCO%id,'load',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%isostasy%load)) then
       call write_log('Creating variable load')
       status = parallel_def_var(NCO%id,'load',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'bedrock deflection from applied load')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable load was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     lon -- longitude
    pos = index(NCO%vars,' lon ')
    status = parallel_inq_varid(NCO%id,'lon',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+3) = '   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%general%lon)) then
       call write_log('Creating variable lon')
       status = parallel_def_var(NCO%id,'lon',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degreeE')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'longitude')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'longitude')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable lon was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     lsurf -- ice lower surface elevation
    pos = index(NCO%vars,' lsurf ')
    status = parallel_inq_varid(NCO%id,'lsurf',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%lsrf)) then
       call write_log('Creating variable lsurf')
       status = parallel_def_var(NCO%id,'lsurf',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice lower surface elevation')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable lsurf was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     marine_connection_mask -- mask for marine_connected cells
    pos = index(NCO%vars,' marine_connection_mask ')
    status = parallel_inq_varid(NCO%id,'marine_connection_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+22) = '                      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%marine_connection_mask)) then
       call write_log('Creating variable marine_connection_mask')
       status = parallel_def_var(NCO%id,'marine_connection_mask',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask for marine_connected cells')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable marine_connection_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     marine_connection_mask_isolated -- mask for marine_connected cells but blocked by grounded cells
    pos = index(NCO%vars,' marine_connection_mask_isolated ')
    status = parallel_inq_varid(NCO%id,'marine_connection_mask_isolated',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+31) = '                               '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%marine_connection_mask_isolated)) then
       call write_log('Creating variable marine_connection_mask_isolated')
       status = parallel_def_var(NCO%id,'marine_connection_mask_isolated',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask for marine_connected cells but blocked by grounded cells')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable marine_connection_mask_isolated was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     overwrite_acab_mask -- cells where acab is overwritten
    pos = index(NCO%vars,' overwrite_acab_mask ')
    status = parallel_inq_varid(NCO%id,'overwrite_acab_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+19) = '                   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%overwrite_acab_mask)) then
       call write_log('Creating variable overwrite_acab_mask')
       status = parallel_def_var(NCO%id,'overwrite_acab_mask',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'cells where acab is overwritten')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_overwrite_acab_mask')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable overwrite_acab_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     powerlaw_c -- spatially varying C for powerlaw sliding, staggered grid
    pos = index(NCO%vars,' powerlaw_c ')
    status = parallel_inq_varid(NCO%id,'powerlaw_c',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%powerlaw_c)) then
       call write_log('Creating variable powerlaw_c')
       status = parallel_def_var(NCO%id,'powerlaw_c',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa (m/yr)**(-1/3)')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'spatially varying C for powerlaw sliding, staggered grid')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable powerlaw_c was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     reference_thck -- reference ice thickness
    pos = index(NCO%vars,' reference_thck ')
    status = parallel_inq_varid(NCO%id,'reference_thck',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%reference_thck)) then
       call write_log('Creating variable reference_thck')
       status = parallel_def_var(NCO%id,'reference_thck',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'reference ice thickness')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable reference_thck was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     relx -- relaxed bedrock topography
    pos = index(NCO%vars,' relx ')
    status = parallel_inq_varid(NCO%id,'relx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%isostasy%relx)) then
       call write_log('Creating variable relx')
       status = parallel_def_var(NCO%id,'relx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'relaxed bedrock topography')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable relx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     resid_u -- u component of residual Ax - b (NOTE: Glissade only)
    pos = index(NCO%vars,' resid_u ')
    status = parallel_inq_varid(NCO%id,'resid_u',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%resid_u)) then
       call write_log('Creating variable resid_u')
       status = parallel_def_var(NCO%id,'resid_u',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_resid))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa/m')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'u component of residual Ax - b (NOTE: Glissade only)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable resid_u was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     resid_v -- v component of residual Ax - b (NOTE: Glissade only)
    pos = index(NCO%vars,' resid_v ')
    status = parallel_inq_varid(NCO%id,'resid_v',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%resid_v)) then
       call write_log('Creating variable resid_v')
       status = parallel_def_var(NCO%id,'resid_v',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_resid))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa/m')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'v component of residual Ax - b (NOTE: Glissade only)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable resid_v was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     rho_ice -- ice density
    pos = index(NCO%vars,' rho_ice ')
    status = parallel_inq_varid(NCO%id,'rho_ice',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(rhoi)) then
       call write_log('Creating variable rho_ice')
       status = parallel_def_var(NCO%id,'rho_ice',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/meter3')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice density')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'rho_ice')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable rho_ice was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     rho_seawater -- seawater density
    pos = index(NCO%vars,' rho_seawater ')
    status = parallel_inq_varid(NCO%id,'rho_seawater',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(rhoo)) then
       call write_log('Creating variable rho_seawater')
       status = parallel_def_var(NCO%id,'rho_seawater',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/meter3')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'seawater density')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'rho_seawater')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable rho_seawater was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     rhs_u -- u component of b in Ax = b
    pos = index(NCO%vars,' rhs_u ')
    status = parallel_inq_varid(NCO%id,'rhs_u',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%rhs_u)) then
       call write_log('Creating variable rhs_u')
       status = parallel_def_var(NCO%id,'rhs_u',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_resid))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa/m')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'u component of b in Ax = b')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable rhs_u was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     rhs_v -- v component of b in Ax = b
    pos = index(NCO%vars,' rhs_v ')
    status = parallel_inq_varid(NCO%id,'rhs_v',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%rhs_v)) then
       call write_log('Creating variable rhs_v')
       status = parallel_def_var(NCO%id,'rhs_v',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_resid))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa/m')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'v component of b in Ax = b')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable rhs_v was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     seconds_per_year -- seconds per year
    pos = index(NCO%vars,' seconds_per_year ')
    status = parallel_inq_varid(NCO%id,'seconds_per_year',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+16) = '                '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(scyr)) then
       call write_log('Creating variable seconds_per_year')
       status = parallel_def_var(NCO%id,'seconds_per_year',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',1.0)
       status = parallel_put_att(NCO%id, varid, 'units', &
            's/yr')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'seconds per year')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'seconds_per_year')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable seconds_per_year was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     sfc_mbal_flux -- surface mass balance flux
    pos = index(NCO%vars,' sfc_mbal_flux ')
    status = parallel_inq_varid(NCO%id,'sfc_mbal_flux',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+13) = '             '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%sfc_mbal_flux)) then
       call write_log('Creating variable sfc_mbal_flux')
       status = parallel_def_var(NCO%id,'sfc_mbal_flux',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/m2/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface mass balance flux')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_flux')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable sfc_mbal_flux was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     sfc_mbal_flux_tavg -- surface mass balance flux (time average)
    pos = index(NCO%vars,' sfc_mbal_flux_tavg ')
    status = parallel_inq_varid(NCO%id,'sfc_mbal_flux_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+18) = '                  '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%sfc_mbal_flux_tavg)) then
       call write_log('Creating variable sfc_mbal_flux_tavg')
       status = parallel_def_var(NCO%id,'sfc_mbal_flux_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/m2/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface mass balance flux (time average)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_flux')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean over years')
       status = parallel_put_att(NCO%id, varid, 'avg_factor', &
            'tavgf')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable sfc_mbal_flux_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     smb -- surface mass balance
    pos = index(NCO%vars,' smb ')
    status = parallel_inq_varid(NCO%id,'smb',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+3) = '   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%smb)) then
       call write_log('Creating variable smb')
       status = parallel_def_var(NCO%id,'smb',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1.0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'mm/year water equivalent')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface mass balance')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable smb was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     smb_3d -- 3D surface mass balance
    pos = index(NCO%vars,' smb_3d ')
    status = parallel_inq_varid(NCO%id,'smb_3d',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%smb_3d)) then
       call write_log('Creating variable smb_3d')
       status = parallel_def_var(NCO%id,'smb_3d',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, nlev_smb_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1.0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'mm/year water equivalent')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            '3D surface mass balance')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_3d')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable smb_3d was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     smb_anomaly -- accumulation, ablation rate anomaly
    pos = index(NCO%vars,' smb_anomaly ')
    status = parallel_inq_varid(NCO%id,'smb_anomaly',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%smb_anomaly)) then
       call write_log('Creating variable smb_anomaly')
       status = parallel_def_var(NCO%id,'smb_anomaly',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1.0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'mm/year water equivalent')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'accumulation, ablation rate anomaly')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_anomaly')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable smb_anomaly was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     smb_gradz -- surface mass balance vertical gradient
    pos = index(NCO%vars,' smb_gradz ')
    status = parallel_inq_varid(NCO%id,'smb_gradz',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%smb_gradz)) then
       call write_log('Creating variable smb_gradz')
       status = parallel_def_var(NCO%id,'smb_gradz',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1.0/thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'mm/year water equivalent per m')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface mass balance vertical gradient')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_vertical_gradient')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable smb_gradz was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     smb_levels -- surface mass balance reference levels
    pos = index(NCO%vars,' smb_levels ')
    status = parallel_inq_varid(NCO%id,'smb_levels',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%smb_levels)) then
       call write_log('Creating variable smb_levels')
       status = parallel_def_var(NCO%id,'smb_levels',get_xtype(outfile,NF90_FLOAT), &
            (/nlev_smb_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'm')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface mass balance reference levels')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_reference_levels')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable smb_levels was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     smb_ref -- surface mass balance at reference elevation
    pos = index(NCO%vars,' smb_ref ')
    status = parallel_inq_varid(NCO%id,'smb_ref',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%smb_ref)) then
       call write_log('Creating variable smb_ref')
       status = parallel_def_var(NCO%id,'smb_ref',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1.0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'mm/year water equivalent')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface mass balance at reference elevation')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_specific_mass_balance_reference')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable smb_ref was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     smb_reference_usrf -- reference upper surface elevation for SMB forcing
    pos = index(NCO%vars,' smb_reference_usrf ')
    status = parallel_inq_varid(NCO%id,'smb_reference_usrf',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+18) = '                  '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%climate%smb_reference_usrf)) then
       call write_log('Creating variable smb_reference_usrf')
       status = parallel_def_var(NCO%id,'smb_reference_usrf',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'm')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'reference upper surface elevation for SMB forcing')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_specific_surface_mass_balance_reference_elevation')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable smb_reference_usrf was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     soft -- bed softness parameter
    pos = index(NCO%vars,' soft ')
    status = parallel_inq_varid(NCO%id,'soft',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%bed_softness)) then
       call write_log('Creating variable soft')
       status = parallel_def_var(NCO%id,'soft',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_btrc))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/pascal/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'bed softness parameter')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable soft was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     stagthk -- staggered ice thickness
    pos = index(NCO%vars,' stagthk ')
    status = parallel_inq_varid(NCO%id,'stagthk',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geomderv%stagthck)) then
       call write_log('Creating variable stagthk')
       status = parallel_def_var(NCO%id,'stagthk',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'staggered ice thickness')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'stag_land_ice_thickness')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable stagthk was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_c -- yield stress for plastic sliding
    pos = index(NCO%vars,' tau_c ')
    status = parallel_inq_varid(NCO%id,'tau_c',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%tau_c)) then
       call write_log('Creating variable tau_c')
       status = parallel_def_var(NCO%id,'tau_c',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1e-3))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kilopascal')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'yield stress for plastic sliding')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_c was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_eff -- effective stress
    pos = index(NCO%vars,' tau_eff ')
    status = parallel_inq_varid(NCO%id,'tau_eff',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%scalar)) then
       call write_log('Creating variable tau_eff')
       status = parallel_def_var(NCO%id,'tau_eff',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'effective stress')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_eff was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_eff_calving -- effective stress for calving
    pos = index(NCO%vars,' tau_eff_calving ')
    status = parallel_inq_varid(NCO%id,'tau_eff_calving',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+15) = '               '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%tau_eff)) then
       call write_log('Creating variable tau_eff_calving')
       status = parallel_def_var(NCO%id,'tau_eff_calving',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'effective stress for calving')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_eff_calving was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_eigen1 -- first eigenvalue of horizontal stress tensor
    pos = index(NCO%vars,' tau_eigen1 ')
    status = parallel_inq_varid(NCO%id,'tau_eigen1',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%tau_eigen1)) then
       call write_log('Creating variable tau_eigen1')
       status = parallel_def_var(NCO%id,'tau_eigen1',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'first eigenvalue of horizontal stress tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_eigen1 was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_eigen2 -- second eigenvalue of horizontal stress tensor
    pos = index(NCO%vars,' tau_eigen2 ')
    status = parallel_inq_varid(NCO%id,'tau_eigen2',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%tau_eigen2)) then
       call write_log('Creating variable tau_eigen2')
       status = parallel_def_var(NCO%id,'tau_eigen2',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'second eigenvalue of horizontal stress tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_eigen2 was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_xx -- xx component of deviatoric stress tensor
    pos = index(NCO%vars,' tau_xx ')
    status = parallel_inq_varid(NCO%id,'tau_xx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%xx)) then
       call write_log('Creating variable tau_xx')
       status = parallel_def_var(NCO%id,'tau_xx',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'xx component of deviatoric stress tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_xx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_xy -- xy component of deviatoric stress tensor
    pos = index(NCO%vars,' tau_xy ')
    status = parallel_inq_varid(NCO%id,'tau_xy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%xy)) then
       call write_log('Creating variable tau_xy')
       status = parallel_def_var(NCO%id,'tau_xy',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'xy component of deviatoric stress tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_xy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_xz -- xz component of deviatoric stress tensor
    pos = index(NCO%vars,' tau_xz ')
    status = parallel_inq_varid(NCO%id,'tau_xz',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%xz)) then
       call write_log('Creating variable tau_xz')
       status = parallel_def_var(NCO%id,'tau_xz',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'xz component of deviatoric stress tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_xz was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_yy -- yy component of deviatoric stress tensor
    pos = index(NCO%vars,' tau_yy ')
    status = parallel_inq_varid(NCO%id,'tau_yy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%yy)) then
       call write_log('Creating variable tau_yy')
       status = parallel_def_var(NCO%id,'tau_yy',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'yy component of deviatoric stress tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_yy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tau_yz -- yz component of deviatoric stress tensor
    pos = index(NCO%vars,' tau_yz ')
    status = parallel_inq_varid(NCO%id,'tau_yz',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+6) = '      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%tau%yz)) then
       call write_log('Creating variable tau_yz')
       status = parallel_def_var(NCO%id,'tau_yz',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'yz component of deviatoric stress tensor')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tau_yz was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     taudx -- driving stress (x-direction comp)
    pos = index(NCO%vars,' taudx ')
    status = parallel_inq_varid(NCO%id,'taudx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%taudx)) then
       call write_log('Creating variable taudx')
       status = parallel_def_var(NCO%id,'taudx',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'driving stress (x-direction comp)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable taudx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     taudy -- driving stress (y-direction comp)
    pos = index(NCO%vars,' taudy ')
    status = parallel_inq_varid(NCO%id,'taudy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%stress%taudy)) then
       call write_log('Creating variable taudy')
       status = parallel_def_var(NCO%id,'taudy',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'driving stress (y-direction comp)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable taudy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tauf -- higher-order basal yield stress
    pos = index(NCO%vars,' tauf ')
    status = parallel_inq_varid(NCO%id,'tauf',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_physics%mintauf)) then
       call write_log('Creating variable tauf')
       status = parallel_def_var(NCO%id,'tauf',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_tau))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'higher-order basal yield stress')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tauf was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     taux -- basal shear stress in x direction (NOTE: Glide only)
    pos = index(NCO%vars,' taux ')
    status = parallel_inq_varid(NCO%id,'taux',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%tau_x)) then
       call write_log('Creating variable taux')
       status = parallel_def_var(NCO%id,'taux',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1e-3*thk0*thk0/len0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kilopascal')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal shear stress in x direction (NOTE: Glide only)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable taux was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tauy -- basal shear stress in y direction
    pos = index(NCO%vars,' tauy ')
    status = parallel_inq_varid(NCO%id,'tauy',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%tau_y)) then
       call write_log('Creating variable tauy')
       status = parallel_def_var(NCO%id,'tauy',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(1e-3*thk0*thk0/len0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kilopascal')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal shear stress in y direction')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tauy was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     temp -- ice temperature
    pos = index(NCO%vars,' temp ')
    status = parallel_inq_varid(NCO%id,'temp',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%temp)) then
       call write_log('Creating variable temp')
       status = parallel_def_var(NCO%id,'temp',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice temperature')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_temperature')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable temp was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tempstag -- ice temperature on staggered vertical levels with boundaries
    pos = index(NCO%vars,' tempstag ')
    status = parallel_inq_varid(NCO%id,'tempstag',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%temp)) then
       call write_log('Creating variable tempstag')
       status = parallel_def_var(NCO%id,'tempstag',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, stagwbndlevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice temperature on staggered vertical levels with boundaries')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_temperature_stag')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tempstag was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     tempunstag -- ice temperature on unstaggered vertical levels
    pos = index(NCO%vars,' tempunstag ')
    status = parallel_inq_varid(NCO%id,'tempunstag',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%tempunstag)) then
       call write_log('Creating variable tempunstag')
       status = parallel_def_var(NCO%id,'tempunstag',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degree_Celsius')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice temperature on unstaggered vertical levels')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_temperature_unstag')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable tempunstag was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     thck_calving_threshold -- thickness threshold for calving ice
    pos = index(NCO%vars,' thck_calving_threshold ')
    status = parallel_inq_varid(NCO%id,'thck_calving_threshold',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+22) = '                      '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%calving%thck_calving_threshold)) then
       call write_log('Creating variable thck_calving_threshold')
       status = parallel_def_var(NCO%id,'thck_calving_threshold',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'thickness threshold for calving ice')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable thck_calving_threshold was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     thermal_forcing -- thermal_forcing
    pos = index(NCO%vars,' thermal_forcing ')
    status = parallel_inq_varid(NCO%id,'thermal_forcing',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+15) = '               '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%ocean_data%thermal_forcing)) then
       call write_log('Creating variable thermal_forcing')
       status = parallel_def_var(NCO%id,'thermal_forcing',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, zocn_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degrees K')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'thermal_forcing')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable thermal_forcing was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     thermal_forcing_lsrf -- thermal_forcing at lower ice surface
    pos = index(NCO%vars,' thermal_forcing_lsrf ')
    status = parallel_inq_varid(NCO%id,'thermal_forcing_lsrf',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+20) = '                    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%ocean_data%thermal_forcing_lsrf)) then
       call write_log('Creating variable thermal_forcing_lsrf')
       status = parallel_def_var(NCO%id,'thermal_forcing_lsrf',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'degrees K')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'thermal_forcing at lower ice surface')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable thermal_forcing_lsrf was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     thk -- ice thickness
    pos = index(NCO%vars,' thk ')
    status = parallel_inq_varid(NCO%id,'thk',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+3) = '   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%thck)) then
       call write_log('Creating variable thk')
       status = parallel_def_var(NCO%id,'thk',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice thickness')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_thickness')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable thk was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     thkmask -- mask
    pos = index(NCO%vars,' thkmask ')
    status = parallel_inq_varid(NCO%id,'thkmask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%thkmask)) then
       call write_log('Creating variable thkmask')
       status = parallel_def_var(NCO%id,'thkmask',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'mask')
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable thkmask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     topg -- bedrock topography
    pos = index(NCO%vars,' topg ')
    status = parallel_inq_varid(NCO%id,'topg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%topg)) then
       call write_log('Creating variable topg')
       status = parallel_def_var(NCO%id,'topg',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'bedrock topography')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'bedrock_altitude')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable topg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     topg_stdev -- standard deviation of bedrock topography
    pos = index(NCO%vars,' topg_stdev ')
    status = parallel_inq_varid(NCO%id,'topg_stdev',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%topg_stdev)) then
       call write_log('Creating variable topg_stdev')
       status = parallel_def_var(NCO%id,'topg_stdev',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'standard deviation of bedrock topography')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'bedrock_altitude_stdev')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable topg_stdev was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     total_bmb_flux -- total basal mass balance flux
    pos = index(NCO%vars,' total_bmb_flux ')
    status = parallel_inq_varid(NCO%id,'total_bmb_flux',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%total_bmb_flux)) then
       call write_log('Creating variable total_bmb_flux')
       status = parallel_def_var(NCO%id,'total_bmb_flux',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'total basal mass balance flux')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable total_bmb_flux was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     total_bmb_flux_tavg -- total basal mass balance flux (time average)
    pos = index(NCO%vars,' total_bmb_flux_tavg ')
    status = parallel_inq_varid(NCO%id,'total_bmb_flux_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+19) = '                   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%total_bmb_flux_tavg)) then
       call write_log('Creating variable total_bmb_flux_tavg')
       status = parallel_def_var(NCO%id,'total_bmb_flux_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'total basal mass balance flux (time average)')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean over years')
       status = parallel_put_att(NCO%id, varid, 'avg_factor', &
            'tavgf')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable total_bmb_flux_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     total_calving_flux -- total calving mass balance flux
    pos = index(NCO%vars,' total_calving_flux ')
    status = parallel_inq_varid(NCO%id,'total_calving_flux',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+18) = '                  '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%total_calving_flux)) then
       call write_log('Creating variable total_calving_flux')
       status = parallel_def_var(NCO%id,'total_calving_flux',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'total calving mass balance flux')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable total_calving_flux was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     total_calving_flux_tavg -- total calving mass balance flux (time average)
    pos = index(NCO%vars,' total_calving_flux_tavg ')
    status = parallel_inq_varid(NCO%id,'total_calving_flux_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+23) = '                       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%total_calving_flux_tavg)) then
       call write_log('Creating variable total_calving_flux_tavg')
       status = parallel_def_var(NCO%id,'total_calving_flux_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'total calving mass balance flux (time average)')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean over years')
       status = parallel_put_att(NCO%id, varid, 'avg_factor', &
            'tavgf')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable total_calving_flux_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     total_gl_flux -- total grounding line flux
    pos = index(NCO%vars,' total_gl_flux ')
    status = parallel_inq_varid(NCO%id,'total_gl_flux',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+13) = '             '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%total_gl_flux)) then
       call write_log('Creating variable total_gl_flux')
       status = parallel_def_var(NCO%id,'total_gl_flux',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'total grounding line flux')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable total_gl_flux was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     total_gl_flux_tavg -- total grounding line flux (time average)
    pos = index(NCO%vars,' total_gl_flux_tavg ')
    status = parallel_inq_varid(NCO%id,'total_gl_flux_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+18) = '                  '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%total_gl_flux_tavg)) then
       call write_log('Creating variable total_gl_flux_tavg')
       status = parallel_def_var(NCO%id,'total_gl_flux_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'total grounding line flux (time average)')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean over years')
       status = parallel_put_att(NCO%id, varid, 'avg_factor', &
            'tavgf')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable total_gl_flux_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     total_smb_flux -- total surface mass balance flux
    pos = index(NCO%vars,' total_smb_flux ')
    status = parallel_inq_varid(NCO%id,'total_smb_flux',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%total_smb_flux)) then
       call write_log('Creating variable total_smb_flux')
       status = parallel_def_var(NCO%id,'total_smb_flux',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'total surface mass balance flux')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable total_smb_flux was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     total_smb_flux_tavg -- total surface mass balance flux (time average)
    pos = index(NCO%vars,' total_smb_flux_tavg ')
    status = parallel_inq_varid(NCO%id,'total_smb_flux_tavg',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+19) = '                   '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%total_smb_flux_tavg)) then
       call write_log('Creating variable total_smb_flux_tavg')
       status = parallel_def_var(NCO%id,'total_smb_flux_tavg',get_xtype(outfile,NF90_FLOAT), &
            (/time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'kg/s')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'total surface mass balance flux (time average)')
       status = parallel_put_att(NCO%id, varid, 'cell_methods', &
            'time: mean over years')
       status = parallel_put_att(NCO%id, varid, 'avg_factor', &
            'tavgf')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable total_smb_flux_tavg was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     ubas -- basal slip velocity in x direction
    pos = index(NCO%vars,' ubas ')
    status = parallel_inq_varid(NCO%id,'ubas',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uvel)) then
       call write_log('Creating variable ubas')
       status = parallel_def_var(NCO%id,'ubas',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal slip velocity in x direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_x_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable ubas was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     uflx -- flux in x direction (NOTE: Glide only)
    pos = index(NCO%vars,' uflx ')
    status = parallel_inq_varid(NCO%id,'uflx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uflx)) then
       call write_log('Creating variable uflx')
       status = parallel_def_var(NCO%id,'uflx',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uflx))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter2/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'flux in x direction (NOTE: Glide only)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable uflx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     unstagbeta -- higher-order bed stress coefficient on the unstaggered grid (NOTE: this will overwrite beta if both are input)
    pos = index(NCO%vars,' unstagbeta ')
    status = parallel_inq_varid(NCO%id,'unstagbeta',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+10) = '          '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%unstagbeta)) then
       call write_log('Creating variable unstagbeta')
       status = parallel_def_var(NCO%id,'unstagbeta',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_beta))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'Pa yr/m')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'higher-order bed stress coefficient on the unstaggered grid (NOTE: this will overwrite beta if both are input)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable unstagbeta was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     usfc -- surface velocity in x direction
    pos = index(NCO%vars,' usfc ')
    status = parallel_inq_varid(NCO%id,'usfc',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uvel)) then
       call write_log('Creating variable usfc')
       status = parallel_def_var(NCO%id,'usfc',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface velocity in x direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_x_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable usfc was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     usfc_obs -- observed surface velocity in x direction
    pos = index(NCO%vars,' usfc_obs ')
    status = parallel_inq_varid(NCO%id,'usfc_obs',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%usfc_obs)) then
       call write_log('Creating variable usfc_obs')
       status = parallel_def_var(NCO%id,'usfc_obs',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'observed surface velocity in x direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_x_velocity_observed')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable usfc_obs was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     usrf_obs -- observed surface elevation
    pos = index(NCO%vars,' usrf_obs ')
    status = parallel_inq_varid(NCO%id,'usrf_obs',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%usrf_obs)) then
       call write_log('Creating variable usrf_obs')
       status = parallel_def_var(NCO%id,'usrf_obs',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'observed surface elevation')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable usrf_obs was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     usurf -- ice upper surface elevation
    pos = index(NCO%vars,' usurf ')
    status = parallel_inq_varid(NCO%id,'usurf',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+5) = '     '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%geometry%usrf)) then
       call write_log('Creating variable usurf')
       status = parallel_def_var(NCO%id,'usurf',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(thk0))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice upper surface elevation')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'surface_altitude')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable usurf was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     uvel -- ice velocity in x direction
    pos = index(NCO%vars,' uvel ')
    status = parallel_inq_varid(NCO%id,'uvel',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uvel)) then
       call write_log('Creating variable uvel')
       status = parallel_def_var(NCO%id,'uvel',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice velocity in x direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_x_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable uvel was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     uvel_2d -- vertically averaged ice velocity in x direction
    pos = index(NCO%vars,' uvel_2d ')
    status = parallel_inq_varid(NCO%id,'uvel_2d',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uvel_2d)) then
       call write_log('Creating variable uvel_2d')
       status = parallel_def_var(NCO%id,'uvel_2d',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertically averaged ice velocity in x direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_x_velocity_2d')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable uvel_2d was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     uvel_2d_extend -- vertically averaged ice velocity in x direction (extended grid)
    pos = index(NCO%vars,' uvel_2d_extend ')
    status = parallel_inq_varid(NCO%id,'uvel_2d_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uvel_2d_extend)) then
       call write_log('Creating variable uvel_2d_extend')
       status = parallel_def_var(NCO%id,'uvel_2d_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertically averaged ice velocity in x direction (extended grid)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_x_velocity_2d')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable uvel_2d_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     uvel_extend -- ice velocity in x direction (extended grid)
    pos = index(NCO%vars,' uvel_extend ')
    status = parallel_inq_varid(NCO%id,'uvel_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uvel_extend)) then
       call write_log('Creating variable uvel_extend')
       status = parallel_def_var(NCO%id,'uvel_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice velocity in x direction (extended grid)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_x_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable uvel_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     uvel_mean -- vertical mean ice velocity in x direction
    pos = index(NCO%vars,' uvel_mean ')
    status = parallel_inq_varid(NCO%id,'uvel_mean',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%uvel_mean)) then
       call write_log('Creating variable uvel_mean')
       status = parallel_def_var(NCO%id,'uvel_mean',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertical mean ice velocity in x direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_vertical_mean_x_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable uvel_mean was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vbas -- basal slip velocity in y direction
    pos = index(NCO%vars,' vbas ')
    status = parallel_inq_varid(NCO%id,'vbas',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vvel)) then
       call write_log('Creating variable vbas')
       status = parallel_def_var(NCO%id,'vbas',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal slip velocity in y direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_y_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vbas was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     velnorm -- Horizontal ice velocity magnitude
    pos = index(NCO%vars,' velnorm ')
    status = parallel_inq_varid(NCO%id,'velnorm',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%velnorm)) then
       call write_log('Creating variable velnorm')
       status = parallel_def_var(NCO%id,'velnorm',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Horizontal ice velocity magnitude')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable velnorm was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     velo_sfc_obs -- observed surface speed
    pos = index(NCO%vars,' velo_sfc_obs ')
    status = parallel_inq_varid(NCO%id,'velo_sfc_obs',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+12) = '            '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%velo_sfc_obs)) then
       call write_log('Creating variable velo_sfc_obs')
       status = parallel_def_var(NCO%id,'velo_sfc_obs',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'observed surface speed')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_speed_observed')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable velo_sfc_obs was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vflx -- flux in x direction (NOTE: Glide only)
    pos = index(NCO%vars,' vflx ')
    status = parallel_inq_varid(NCO%id,'vflx',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vflx)) then
       call write_log('Creating variable vflx')
       status = parallel_def_var(NCO%id,'vflx',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uflx))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter2/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'flux in x direction (NOTE: Glide only)')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vflx was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vsfc -- surface velocity in y direction
    pos = index(NCO%vars,' vsfc ')
    status = parallel_inq_varid(NCO%id,'vsfc',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vvel)) then
       call write_log('Creating variable vsfc')
       status = parallel_def_var(NCO%id,'vsfc',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface velocity in y direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_y_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vsfc was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vsfc_obs -- observed surface velocity in y direction
    pos = index(NCO%vars,' vsfc_obs ')
    status = parallel_inq_varid(NCO%id,'vsfc_obs',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+8) = '        '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vsfc_obs)) then
       call write_log('Creating variable vsfc_obs')
       status = parallel_def_var(NCO%id,'vsfc_obs',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'observed surface velocity in y direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_y_velocity_observed')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vsfc_obs was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vvel -- ice velocity in y direction
    pos = index(NCO%vars,' vvel ')
    status = parallel_inq_varid(NCO%id,'vvel',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vvel)) then
       call write_log('Creating variable vvel')
       status = parallel_def_var(NCO%id,'vvel',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice velocity in y direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_y_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vvel was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vvel_2d -- vertically averaged ice velocity in y direction
    pos = index(NCO%vars,' vvel_2d ')
    status = parallel_inq_varid(NCO%id,'vvel_2d',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+7) = '       '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vvel_2d)) then
       call write_log('Creating variable vvel_2d')
       status = parallel_def_var(NCO%id,'vvel_2d',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertically averaged ice velocity in y direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_y_velocity_2d')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vvel_2d was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vvel_2d_extend -- vertically averaged ice velocity in y direction (extended grid)
    pos = index(NCO%vars,' vvel_2d_extend ')
    status = parallel_inq_varid(NCO%id,'vvel_2d_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+14) = '              '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vvel_2d_extend)) then
       call write_log('Creating variable vvel_2d_extend')
       status = parallel_def_var(NCO%id,'vvel_2d_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertically averaged ice velocity in y direction (extended grid)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_y_velocity_2d')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vvel_2d_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vvel_extend -- ice velocity in y direction (extended grid)
    pos = index(NCO%vars,' vvel_extend ')
    status = parallel_inq_varid(NCO%id,'vvel_extend',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+11) = '           '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vvel_extend)) then
       call write_log('Creating variable vvel_extend')
       status = parallel_def_var(NCO%id,'vvel_extend',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'ice velocity in y direction (extended grid)')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_y_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vvel_extend was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     vvel_mean -- vertical mean ice velocity in y direction
    pos = index(NCO%vars,' vvel_mean ')
    status = parallel_inq_varid(NCO%id,'vvel_mean',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%vvel_mean)) then
       call write_log('Creating variable vvel_mean')
       status = parallel_def_var(NCO%id,'vvel_mean',get_xtype(outfile,NF90_FLOAT), &
            (/x0_dimid, y0_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertical mean ice velocity in y direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_vertical_mean_y_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable vvel_mean was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     warm_ocean_mask -- warm ocean mask
    pos = index(NCO%vars,' warm_ocean_mask ')
    status = parallel_inq_varid(NCO%id,'warm_ocean_mask',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+15) = '               '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%basal_melt%warm_ocean_mask)) then
       call write_log('Creating variable warm_ocean_mask')
       status = parallel_def_var(NCO%id,'warm_ocean_mask',get_xtype(outfile,NF90_INT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            '1')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'warm ocean mask')
       status = parallel_put_att(NCO%id, varid, '_FillValue', -2147483647)
       status = parallel_put_att(NCO%id, varid, 'missing_value', -2147483647)
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable warm_ocean_mask was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     waterfrac -- internal water fraction
    pos = index(NCO%vars,' waterfrac ')
    status = parallel_inq_varid(NCO%id,'waterfrac',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+9) = '         '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%temper%waterfrac)) then
       call write_log('Creating variable waterfrac')
       status = parallel_def_var(NCO%id,'waterfrac',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, staglevel_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'units', &
            'unitless [0,1]')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'internal water fraction')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable waterfrac was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     wbas -- basal velocity in z direction
    pos = index(NCO%vars,' wbas ')
    status = parallel_inq_varid(NCO%id,'wbas',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%wvel)) then
       call write_log('Creating variable wbas')
       status = parallel_def_var(NCO%id,'wbas',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'basal velocity in z direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_basal_z_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable wbas was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     wgrd -- Vertical grid velocity
    pos = index(NCO%vars,' wgrd ')
    status = parallel_inq_varid(NCO%id,'wgrd',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%wgrd)) then
       call write_log('Creating variable wgrd')
       status = parallel_def_var(NCO%id,'wgrd',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_wvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'Vertical grid velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable wgrd was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     wsfc -- surface velocity in z direction
    pos = index(NCO%vars,' wsfc ')
    status = parallel_inq_varid(NCO%id,'wsfc',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%wvel)) then
       call write_log('Creating variable wsfc')
       status = parallel_def_var(NCO%id,'wsfc',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_uvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'surface velocity in z direction')
       status = parallel_put_att(NCO%id, varid, 'standard_name', &
            'land_ice_surface_z_velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable wsfc was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

    !     wvel -- vertical ice velocity
    pos = index(NCO%vars,' wvel ')
    status = parallel_inq_varid(NCO%id,'wvel',varid)
    if (pos.ne.0) then
      NCO%vars(pos+1:pos+4) = '    '
    end if
    if (pos.ne.0 .and. status.eq.nf90_enotvar) then
    if (is_enabled(data%velocity%wvel)) then
       call write_log('Creating variable wvel')
       status = parallel_def_var(NCO%id,'wvel',get_xtype(outfile,NF90_FLOAT), &
            (/x1_dimid, y1_dimid, level_dimid, time_dimid/),varid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, varid, 'scale_factor',(scale_wvel))
       status = parallel_put_att(NCO%id, varid, 'units', &
            'meter/year')
       status = parallel_put_att(NCO%id, varid, 'long_name', &
            'vertical ice velocity')
       if (get_xtype(outfile,NF90_FLOAT) == NF90_DOUBLE) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0d20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0d20)
       elseif (get_xtype(outfile,NF90_FLOAT) == NF90_FLOAT) then
          status = parallel_put_att(NCO%id, varid, '_FillValue', 1.0e20)
          status = parallel_put_att(NCO%id, varid, 'missing_value', 1.0e20)
       endif
       if (glimmap_allocated(model%projection)) then
          status = parallel_put_att(NCO%id, varid, 'grid_mapping',glimmer_nc_mapvarname)
       end if
     else
     call write_log('Variable wvel was specified for output but it is &
          &inappropriate for your config settings.  It will be excluded from the output.', GM_WARNING)
     end if
     end if

  end subroutine glide_io_create

  subroutine glide_io_write(outfile,data)

    use cism_parallel, only: parallel_type, parallel_inq_varid, distributed_put_var, parallel_put_var
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: data
    ! the model instance

    ! local variables
    real(dp) :: tavgf
    integer status, varid
    integer up
    type(parallel_type) :: parallel

    parallel = data%parallel
     
    tavgf = outfile%total_time
    if (tavgf.ne.0.d0) then
       tavgf = 1.d0/tavgf
    end if

    ! write variables
    status = parallel_inq_varid(NCO%id,'S_ambient',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%plume%S_ambient(1:data%general%ewn,1:data%general%nsn), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'T_ambient',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%plume%T_ambient(1:data%general%ewn,1:data%general%nsn), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'acab',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%acab, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'acab_3d',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%acab_3d(up,:,:), parallel, (/1,1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'acab_anomaly',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%acab_anomaly, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'acab_applied',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%acab_applied, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'acab_applied_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            (tavgf)*(data%climate%acab_applied_tavg), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'acab_corrected',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%acab_corrected, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'acab_gradz',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%acab_gradz, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'acab_ref',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%acab_ref, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'adv_cfl_dt',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%numerics%adv_cfl_dt, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'area_factor',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%projection%stere%area_factor, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'artm',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%artm, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'artm_3d',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%artm_3d(up,:,:), parallel, (/1,1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'artm_anomaly',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%artm_anomaly, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'artm_gradz',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%artm_gradz, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'artm_ref',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%artm_ref, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'basal_mbal_flux',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%basal_mbal_flux, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'basal_mbal_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            (tavgf)*(data%geometry%basal_mbal_flux_tavg), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'basin_multiplier_array',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_melt%basin_multiplier_array, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'basin_number',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%ocean_data%basin_number, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'beta',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%beta, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'beta_internal',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%beta_internal, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bfricflx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%bfricflx, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bheatflx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%bheatflx, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bmlt',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_melt%bmlt, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bmlt_applied',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_melt%bmlt_applied, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bmlt_applied_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            (tavgf)*(data%basal_melt%bmlt_applied_tavg), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bmlt_float',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_melt%bmlt_float, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bmlt_float_anomaly',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_melt%bmlt_float_anomaly, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bmlt_float_external',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_melt%bmlt_float_external, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bmlt_ground',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_melt%bmlt_ground, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bpmp',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%bpmp, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btemp',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%btemp, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btemp_float',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%btemp_float, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btemp_ground',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%btemp_ground, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btract',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%btract(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btractx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%btractx(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btractx_extend',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%btractx_extend(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btracty',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%btracty(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btracty_extend',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%btracty_extend(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'btrc',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%btrc, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bwat',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_hydro%bwat, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bwat_diag',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_hydro%bwat_diag, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'bwatflx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_hydro%bwatflx, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'c_flux_array',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_hydro%c_flux_array, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'c_space_factor',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%c_space_factor, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'calving_flux',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%calving_flux, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'calving_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            (tavgf)*(data%geometry%calving_flux_tavg), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'calving_lateral',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%lateral_rate, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'calving_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%calving_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'calving_rate',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%calving_rate, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'calving_rate_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            (tavgf)*(data%calving%calving_rate_tavg), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'calving_thck',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%calving_thck, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'cell_area',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%cell_area, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'coulomb_c',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%coulomb_c, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'coulomb_c_relax',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%coulomb_c_relax, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'damage',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%calving%damage(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'deltaT_ocn',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%ocean_data%deltaT_ocn, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'diff_cfl_dt',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%numerics%diff_cfl_dt, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'diffu',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%diffu, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'dissip',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%dissip(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'dissipstag',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%dissip(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'divu',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%divu, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'dthck_dt',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%dthck_dt, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'dthck_dt_obs',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%dthck_dt_obs, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'dthck_dt_obs_basin',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%dthck_dt_obs_basin, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'dthckdtm',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geomderv%dthckdtm, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'dusrfdtm',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geomderv%dusrfdtm, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'effecpress',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%effecpress, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'efvs',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%efvs(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'enthalpy',varid)
    if (status .eq. nf90_noerr) then
       do up=0,NCO%nstagwbndlevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%enthalpy(up,:,:), parallel, (/1,1,up+1,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'eps_eff',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%strain_rate%scalar(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'eps_eigen1',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%eps_eigen1, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'eps_eigen2',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%eps_eigen2, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'eps_xx',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%strain_rate%xx(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'eps_xy',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%strain_rate%xy(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'eps_xz',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%strain_rate%xz(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'eps_yy',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%strain_rate%yy(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'eps_yz',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%strain_rate%yz(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'eus',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%climate%eus, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'f_effecpress_bwat',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%f_effecpress_bwat, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'f_effecpress_bwat_target',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%f_effecpress_bwat_target, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'f_effecpress_bwatflx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%f_effecpress_bwatflx, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'f_effecpress_ocean_p',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%f_effecpress_ocean_p, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'f_flotation',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%f_flotation, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'f_ground',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%f_ground, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'f_ground_cell',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%f_ground_cell, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'f_ground_obs',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%f_ground_obs, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ff_invert_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%inversion%ff_invert_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'floating_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%floating_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'floating_thck_target',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%inversion%floating_thck_target, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'flow_enhancement_factor',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%temper%flow_enhancement_factor, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'flwa',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%flwa(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'flwastag',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%flwa(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'gl_flux',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%gl_flux, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'gl_flux_east',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%gl_flux_east, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'gl_flux_north',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%gl_flux_north, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'gl_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            (tavgf)*(data%geometry%gl_flux_tavg), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'gravity',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            grav, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'grounded_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%grounded_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'head',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_hydro%head, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'iarea',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%geometry%iarea, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'iareaf',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%geometry%iareaf, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'iareag',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%geometry%iareag, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ice_age',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%geometry%ice_age(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'ice_cap_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%ice_cap_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ice_domain_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%general%ice_domain_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ice_fraction_retreat_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%ice_fraction_retreat_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ice_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%ice_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ice_mask_stag',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%ice_mask_stag, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ice_sheet_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%ice_sheet_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ice_specific_heat',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            shci, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ice_thermal_conductivity',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            coni, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'imass',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%geometry%imass, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'imass_above_flotation',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%geometry%imass_above_flotation, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ivol',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%geometry%ivol, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'kinbcmask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%kinbcmask(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'lat',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%general%lat, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'litho_temp',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%lithot%temp, parallel, (/1,1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'load',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%isostasy%load, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'lon',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%general%lon, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'lsurf',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%lsrf, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'marine_connection_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%marine_connection_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'marine_connection_mask_isolated',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%marine_connection_mask_isolated, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'overwrite_acab_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%overwrite_acab_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'powerlaw_c',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%powerlaw_c, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'reference_thck',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%reference_thck, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'relx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%isostasy%relx, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'resid_u',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%resid_u(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'resid_v',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%resid_v(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'rho_ice',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            rhoi, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'rho_seawater',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            rhoo, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'rhs_u',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%rhs_u(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'rhs_v',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%rhs_v(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'seconds_per_year',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            scyr, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'sfc_mbal_flux',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%sfc_mbal_flux, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'sfc_mbal_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            (tavgf)*(data%geometry%sfc_mbal_flux_tavg), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'smb',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%smb, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'smb_3d',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%smb_3d(up,:,:), parallel, (/1,1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'smb_anomaly',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%smb_anomaly, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'smb_gradz',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%smb_gradz, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'smb_levels',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%smb_levels, parallel, (/1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'smb_ref',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%smb_ref, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'smb_reference_usrf',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%climate%smb_reference_usrf, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'soft',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%bed_softness, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'stagthk',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geomderv%stagthck, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'tau_c',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%tau_c, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'tau_eff',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%scalar(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tau_eff_calving',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%tau_eff, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'tau_eigen1',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%tau_eigen1, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'tau_eigen2',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%tau_eigen2, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'tau_xx',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%xx(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tau_xy',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%xy(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tau_xz',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%xz(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tau_yy',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%yy(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tau_yz',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%stress%tau%yz(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'taudx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%taudx(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'taudy',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%stress%taudy(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'tauf',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_physics%mintauf, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'taux',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%tau_x, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'tauy',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%tau_y, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'temp',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%temp(up,1:data%general%ewn,1:data%general%nsn), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tempstag',varid)
    if (status .eq. nf90_noerr) then
       do up=0,NCO%nstagwbndlevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%temp(up,1:data%general%ewn,1:data%general%nsn), parallel, (/1,1,up+1,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'tempunstag',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%tempunstag(up,1:data%general%ewn,1:data%general%nsn), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'thck_calving_threshold',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%calving%thck_calving_threshold, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'thermal_forcing',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nzocn
          status = distributed_put_var(NCO%id, varid, &
               data%ocean_data%thermal_forcing(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'thermal_forcing_lsrf',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%ocean_data%thermal_forcing_lsrf(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'thk',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%thck, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'thkmask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%thkmask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'topg',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%topg, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'topg_stdev',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%topg_stdev, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'total_bmb_flux',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%geometry%total_bmb_flux, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'total_bmb_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            (tavgf)*(data%geometry%total_bmb_flux_tavg), (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'total_calving_flux',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%geometry%total_calving_flux, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'total_calving_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            (tavgf)*(data%geometry%total_calving_flux_tavg), (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'total_gl_flux',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%geometry%total_gl_flux, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'total_gl_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            (tavgf)*(data%geometry%total_gl_flux_tavg), (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'total_smb_flux',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            data%geometry%total_smb_flux, (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'total_smb_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       status = parallel_put_var(NCO%id, varid, &
            (tavgf)*(data%geometry%total_smb_flux_tavg), (/outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'ubas',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%uvel(data%general%upn,1:data%general%ewn-1,1:data%general%nsn-1), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'uflx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%uflx, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'unstagbeta',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%unstagbeta, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'usfc',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%uvel(1,1:data%general%ewn-1,1:data%general%nsn-1), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'usfc_obs',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%usfc_obs, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'usrf_obs',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%usrf_obs, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'usurf',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%geometry%usrf, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'uvel',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%uvel(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'uvel_2d',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%uvel_2d(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'uvel_2d_extend',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%uvel_2d_extend(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'uvel_extend',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%uvel_extend(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'uvel_mean',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%uvel_mean, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'vbas',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%vvel(data%general%upn,1:data%general%ewn-1,1:data%general%nsn-1), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'velnorm',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%velnorm(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'velo_sfc_obs',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%velo_sfc_obs, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'vflx',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%vflx, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'vsfc',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%vvel(1,1:data%general%ewn-1,1:data%general%nsn-1), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'vsfc_obs',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%vsfc_obs, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'vvel',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%vvel(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'vvel_2d',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%vvel_2d(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'vvel_2d_extend',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%vvel_2d_extend(:,:), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'vvel_extend',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%vvel_extend(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'vvel_mean',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%vvel_mean, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'warm_ocean_mask',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%basal_melt%warm_ocean_mask, parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'waterfrac',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nstaglevel
          status = distributed_put_var(NCO%id, varid, &
               data%temper%waterfrac(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'wbas',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%wvel(data%general%upn,1:data%general%ewn,1:data%general%nsn), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'wgrd',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%wgrd(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if

    status = parallel_inq_varid(NCO%id,'wsfc',varid)
    if (status .eq. nf90_noerr) then
       status = distributed_put_var(NCO%id, varid, &
            data%velocity%wvel(1,1:data%general%ewn,1:data%general%nsn), parallel, (/1,1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if

    status = parallel_inq_varid(NCO%id,'wvel',varid)
    if (status .eq. nf90_noerr) then
       do up=1,NCO%nlevel
          status = distributed_put_var(NCO%id, varid, &
               data%velocity%wvel(up,:,:), parallel, (/1,1,up,outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end do
    end if


  end subroutine glide_io_write


  subroutine glide_add_to_restart_variable_list(vars_to_add)
    ! This subroutine adds variables to the list of variables needed for a restart.
    ! It is a public subroutine that allows other parts of the model to modify the list, 
    ! which is a module level variable.   MJH 1/17/2013

    use glimmer_log
    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    character(len=*), intent (in) :: vars_to_add  ! list of variable(s) to be added to the list of restart variables 
    !character(*), intent (inout) :: restart_variable_list  ! list of variables needed to perform an exact restart - module variable

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    ! Add the variables to the list so long as they don't make the list too long.
    if ( (len_trim(restart_variable_list) + 1 + len_trim(vars_to_add)) > len(restart_variable_list)) then
       call write_log('Adding restart variables has made the restart variable list too long.',GM_FATAL)
    else
       restart_variable_list = trim(adjustl(restart_variable_list)) // ' ' // trim(vars_to_add)
       !call write_log('Adding to glide restart variable list: ' // trim(vars_to_add) )
    endif

  end subroutine glide_add_to_restart_variable_list


  ! Functions for the interface 'is_enabled'.  These are needed by the auto-generated code in glide_io_create
  !   to determine if a variable is 'turned on', and should be written.

  function is_enabled_0dint(var)
    integer, intent(in) :: var
    logical :: is_enabled_0dint
    is_enabled_0dint = .true.  ! scalars are always enabled
    return
  end function is_enabled_0dint

  function is_enabled_1dint(var)
    integer, dimension(:), pointer, intent(in) :: var
    logical :: is_enabled_1dint
    if (associated(var)) then
      is_enabled_1dint = .true.
    else
      is_enabled_1dint = .false.
    endif
    return
  end function is_enabled_1dint

  function is_enabled_2dint(var)
    integer, dimension(:,:), pointer, intent(in) :: var
    logical :: is_enabled_2dint
    if (associated(var)) then
      is_enabled_2dint = .true.
    else
      is_enabled_2dint = .false.
    endif
    return
  end function is_enabled_2dint

  function is_enabled_0dreal(var)
    real(dp), intent(in) :: var
    logical :: is_enabled_0dreal
    is_enabled_0dreal = .true.  ! scalars are always enabled
    return
  end function is_enabled_0dreal

  function is_enabled_1dreal(var)
    real(dp), dimension(:), pointer, intent(in) :: var
    logical :: is_enabled_1dreal
    if (associated(var)) then
      is_enabled_1dreal = .true.
    else
      is_enabled_1dreal = .false.
    endif
    return
  end function is_enabled_1dreal

  function is_enabled_2dreal(var)
    real(dp), dimension(:,:), pointer, intent(in) :: var
    logical :: is_enabled_2dreal
    if (associated(var)) then
      is_enabled_2dreal = .true.
    else
      is_enabled_2dreal = .false.
    endif
    return
  end function is_enabled_2dreal

  function is_enabled_3dreal(var)
    real(dp), dimension(:,:,:), pointer, intent(in) :: var
    logical :: is_enabled_3dreal
    if (associated(var)) then
      is_enabled_3dreal = .true.
    else
      is_enabled_3dreal = .false.
    endif
    return
  end function is_enabled_3dreal


  !*****************************************************************************
  ! netCDF input
  !*****************************************************************************  
  subroutine glide_io_readall(data, model, filetype)
    ! read from netCDF file
    use glimmer_ncio
    implicit none
    type(glide_global_type) :: data
    type(glide_global_type) :: model
    integer, intent(in), optional :: filetype  ! 0 for input, 1 for forcing; defaults to input

    ! local variables
    type(glimmer_nc_input), pointer :: ic
    integer :: filetype_local

    if (present(filetype)) then
      filetype_local = filetype
    else
      filetype_local = 0 ! default to input type
    end if

    if (filetype_local == 0) then
      ic=>model%funits%in_first
    else
      ic=>model%funits%frc_first
    endif
    do while(associated(ic))
       call glimmer_nc_checkread(ic,model)
       if (ic%nc%just_processed) then
          call glide_io_read(ic,data)
       end if
       ic=>ic%next
    end do

  end subroutine glide_io_readall


  subroutine glide_read_forcing(data, model)

    ! Read data from forcing files
    use glimmer_log
    use cism_parallel, only: main_task

    implicit none
    type(glide_global_type) :: data
    type(glide_global_type), intent(inout) :: model

    ! Locals
    type(glimmer_nc_input), pointer :: ic
    integer :: t, t_prev
    real(dp) :: current_forcing_time   ! current time with reference to the forcing file
    real(dp) :: eps                    ! a tolerance to use for stepwise constant forcing
    logical, parameter :: verbose_read_forcing = .false.

    ! Make eps a fraction of the time step.
    eps = model%numerics%tinc * 1.0d-4

    ! read forcing files
    ic=>model%funits%frc_first
    do while(associated(ic))

!       if (main_task .and. verbose_read_forcing) print *, 'possible forcing times', ic%times

       ic%nc%just_processed = .true. ! until we find an acceptable time, set this to true which will prevent the file from being read.

       ! Compute the current forcing time.
       ! This is the current model time, plus any offset to be consistent with the time in the forcing file,
       !  plus a small number to allow for roundoff error.
       current_forcing_time = model%numerics%time + ic%time_offset + eps

       ! If cycling repeatedly through a subset of the forcing data, make a further correction:
       ! compute the current time relative to time_start_cycle.
       if (ic%nyear_cycle > 0 .and. current_forcing_time > ic%time_start_cycle) then
          current_forcing_time = ic%time_start_cycle &
               + mod(current_forcing_time - ic%time_start_cycle, real(ic%nyear_cycle,dp))
       endif

       if (main_task .and. verbose_read_forcing) then
          print*, 'In glide_read_forcing, model time + eps =', model%numerics%time + eps
          print*, 'Forcing file nt, time_offset =', ic%nt, ic%time_offset
          print*, 'time_start_cycle, nyear_cycle:', ic%time_start_cycle, ic%nyear_cycle
          print*, 'current forcing time =', current_forcing_time
       endif

       ! Find the time index associated with the previous model time step
       t_prev = 0
       do t = ic%nt, 1, -1  ! look through the time array backwards
          if (ic%times(t) <= current_forcing_time - model%numerics%tinc) then
             t_prev = t
             if (main_task .and. verbose_read_forcing) print*, 'Previous time index =', t_prev
             exit
          end if
       enddo

       ! Find the current time in the file
       do t = ic%nt, 1, -1  ! look through the time array backwards
          if ( ic%times(t) <= current_forcing_time) then
             ! use the largest time that is smaller or equal to the current time (stepwise forcing)
             if (main_task .and. verbose_read_forcing) &
                  print*, 'Largest time less than current forcing time: t, times(t):', t, ic%times(t)

             ! If this time index (t) is larger than the previous index (t_prev), then read a new time slice.
             ! Otherwise, we already have the current slice, and there is nothing new to read.
             if (t > t_prev) then
                ! Set the desired time to be read
                ic%current_time = t
                ic%nc%just_processed = .false.  ! set this to false so file will be read.
                if (main_task .and. verbose_read_forcing) print*, 'Read new forcing slice: t, times(t) =', t, ic%times(t)
             endif ! t > t_prev

             exit  ! once we find the time, exit the loop
          end if   ! ic%times(t) <= model%numerics%time + eps

       end do  ! if we get to end of loop without exiting, then this file will not be read at this time

       ! move on to the next forcing file
       ic=>ic%next
    end do

    ! Now that we've updated metadata for each forcing file, actually perform the read.
    ! This call will only read forcing files where just_processed=.false.
    call glide_io_readall(data, model, filetype=1)

  end subroutine glide_read_forcing


!------------------------------------------------------------------------------


  subroutine glide_io_read(infile,data)

    ! read variables from a netCDF file
    use cism_parallel, only: parallel_type, &
         parallel_inq_varid, parallel_get_att, distributed_get_var, parallel_get_var
    use glimmer_log
    implicit none
    type(glimmer_nc_input), pointer :: infile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: data
    ! the model instance

    ! local variables
    integer status,varid
    integer up
    real(dp) :: scaling_factor
    type(parallel_type) :: parallel

    parallel = data%parallel

    ! read variables
    status = parallel_inq_varid(NCI%id,'x1',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%general%x1)) then
       call write_log('  Loading x1')
       status = distributed_get_var(NCI%id, varid, &
            data%general%x1, parallel, (/1/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling x1",GM_DIAGNOSTIC)
          data%general%x1 = &
               data%general%x1*scaling_factor
       end if
       ! Also read this variable into a global array
       status = parallel_get_var(NCI%id, varid, &
            data%general%x1_global)
       call nc_errorhandle(__FILE__,__LINE__,status)
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling x1_global",GM_DIAGNOSTIC)
          data%general%x1_global = &
               data%general%x1_global*scaling_factor
       end if
    else
    call write_log('Variable x1 was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'y1',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%general%y1)) then
       call write_log('  Loading y1')
       status = distributed_get_var(NCI%id, varid, &
            data%general%y1, parallel, (/1/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling y1",GM_DIAGNOSTIC)
          data%general%y1 = &
               data%general%y1*scaling_factor
       end if
       ! Also read this variable into a global array
       status = parallel_get_var(NCI%id, varid, &
            data%general%y1_global)
       call nc_errorhandle(__FILE__,__LINE__,status)
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling y1_global",GM_DIAGNOSTIC)
          data%general%y1_global = &
               data%general%y1_global*scaling_factor
       end if
    else
    call write_log('Variable y1 was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'acab',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%acab)) then
       call write_log('  Loading acab')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%acab, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_acab)
       else
          scaling_factor = scaling_factor/(scale_acab)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling acab",GM_DIAGNOSTIC)
          data%climate%acab = &
               data%climate%acab*scaling_factor
       end if
    else
    call write_log('Variable acab was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'acab_3d',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%acab_3d)) then
       call write_log('  Loading acab_3d')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%acab_3d(up,:,:), parallel, (/1,1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_acab)
       else
          scaling_factor = scaling_factor/(scale_acab)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling acab_3d",GM_DIAGNOSTIC)
          data%climate%acab_3d(up,:,:) = &
               data%climate%acab_3d(up,:,:)*scaling_factor
       end if
    else
    call write_log('Variable acab_3d was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'acab_anomaly',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%acab_anomaly)) then
       call write_log('  Loading acab_anomaly')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%acab_anomaly, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_acab)
       else
          scaling_factor = scaling_factor/(scale_acab)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling acab_anomaly",GM_DIAGNOSTIC)
          data%climate%acab_anomaly = &
               data%climate%acab_anomaly*scaling_factor
       end if
    else
    call write_log('Variable acab_anomaly was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'acab_gradz',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%acab_gradz)) then
       call write_log('  Loading acab_gradz')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%acab_gradz, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_acab/thk0)
       else
          scaling_factor = scaling_factor/(scale_acab/thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling acab_gradz",GM_DIAGNOSTIC)
          data%climate%acab_gradz = &
               data%climate%acab_gradz*scaling_factor
       end if
    else
    call write_log('Variable acab_gradz was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'acab_ref',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%acab_ref)) then
       call write_log('  Loading acab_ref')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%acab_ref, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_acab)
       else
          scaling_factor = scaling_factor/(scale_acab)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling acab_ref",GM_DIAGNOSTIC)
          data%climate%acab_ref = &
               data%climate%acab_ref*scaling_factor
       end if
    else
    call write_log('Variable acab_ref was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'artm',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%artm)) then
       call write_log('  Loading artm')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%artm, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling artm",GM_DIAGNOSTIC)
          data%climate%artm = &
               data%climate%artm*scaling_factor
       end if
    else
    call write_log('Variable artm was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'artm_3d',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%artm_3d)) then
       call write_log('  Loading artm_3d')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%artm_3d(up,:,:), parallel, (/1,1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling artm_3d",GM_DIAGNOSTIC)
          data%climate%artm_3d(up,:,:) = &
               data%climate%artm_3d(up,:,:)*scaling_factor
       end if
    else
    call write_log('Variable artm_3d was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'artm_anomaly',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%artm_anomaly)) then
       call write_log('  Loading artm_anomaly')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%artm_anomaly, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling artm_anomaly",GM_DIAGNOSTIC)
          data%climate%artm_anomaly = &
               data%climate%artm_anomaly*scaling_factor
       end if
    else
    call write_log('Variable artm_anomaly was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'artm_gradz',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%artm_gradz)) then
       call write_log('  Loading artm_gradz')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%artm_gradz, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(1./thk0)
       else
          scaling_factor = scaling_factor/(1./thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling artm_gradz",GM_DIAGNOSTIC)
          data%climate%artm_gradz = &
               data%climate%artm_gradz*scaling_factor
       end if
    else
    call write_log('Variable artm_gradz was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'artm_ref',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%artm_ref)) then
       call write_log('  Loading artm_ref')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%artm_ref, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling artm_ref",GM_DIAGNOSTIC)
          data%climate%artm_ref = &
               data%climate%artm_ref*scaling_factor
       end if
    else
    call write_log('Variable artm_ref was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'basin_multiplier_array',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_melt%basin_multiplier_array)) then
       call write_log('  Loading basin_multiplier_array')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_melt%basin_multiplier_array, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling basin_multiplier_array",GM_DIAGNOSTIC)
          data%basal_melt%basin_multiplier_array = &
               data%basal_melt%basin_multiplier_array*scaling_factor
       end if
    else
    call write_log('Variable basin_multiplier_array was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'basin_number',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%ocean_data%basin_number)) then
       call write_log('  Loading basin_number')
       status = distributed_get_var(NCI%id, varid, &
            data%ocean_data%basin_number, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling basin_number",GM_DIAGNOSTIC)
          data%ocean_data%basin_number = &
               data%ocean_data%basin_number*scaling_factor
       end if
    else
    call write_log('Variable basin_number was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'beta',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%beta)) then
       call write_log('  Loading beta')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%beta, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_beta)
       else
          scaling_factor = scaling_factor/(scale_beta)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling beta",GM_DIAGNOSTIC)
          data%velocity%beta = &
               data%velocity%beta*scaling_factor
       end if
    else
    call write_log('Variable beta was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bfricflx',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%bfricflx)) then
       call write_log('  Loading bfricflx')
       status = distributed_get_var(NCI%id, varid, &
            data%temper%bfricflx, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bfricflx",GM_DIAGNOSTIC)
          data%temper%bfricflx = &
               data%temper%bfricflx*scaling_factor
       end if
    else
    call write_log('Variable bfricflx was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bheatflx',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%bheatflx)) then
       call write_log('  Loading bheatflx')
       status = distributed_get_var(NCI%id, varid, &
            data%temper%bheatflx, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bheatflx",GM_DIAGNOSTIC)
          data%temper%bheatflx = &
               data%temper%bheatflx*scaling_factor
       end if
    else
    call write_log('Variable bheatflx was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bmlt_float_anomaly',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_melt%bmlt_float_anomaly)) then
       call write_log('  Loading bmlt_float_anomaly')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_melt%bmlt_float_anomaly, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_acab)
       else
          scaling_factor = scaling_factor/(scale_acab)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bmlt_float_anomaly",GM_DIAGNOSTIC)
          data%basal_melt%bmlt_float_anomaly = &
               data%basal_melt%bmlt_float_anomaly*scaling_factor
       end if
    else
    call write_log('Variable bmlt_float_anomaly was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bmlt_float_external',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_melt%bmlt_float_external)) then
       call write_log('  Loading bmlt_float_external')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_melt%bmlt_float_external, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_acab)
       else
          scaling_factor = scaling_factor/(scale_acab)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bmlt_float_external",GM_DIAGNOSTIC)
          data%basal_melt%bmlt_float_external = &
               data%basal_melt%bmlt_float_external*scaling_factor
       end if
    else
    call write_log('Variable bmlt_float_external was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'btract',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%btract)) then
       call write_log('  Loading btract')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%btract(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling btract",GM_DIAGNOSTIC)
          data%stress%btract(:,:) = &
               data%stress%btract(:,:)*scaling_factor
       end if
    else
    call write_log('Variable btract was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'btractx',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%btractx)) then
       call write_log('  Loading btractx')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%btractx(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling btractx",GM_DIAGNOSTIC)
          data%stress%btractx(:,:) = &
               data%stress%btractx(:,:)*scaling_factor
       end if
    else
    call write_log('Variable btractx was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'btractx_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%btractx_extend)) then
       call write_log('  Loading btractx_extend')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%btractx_extend(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling btractx_extend",GM_DIAGNOSTIC)
          data%stress%btractx_extend(:,:) = &
               data%stress%btractx_extend(:,:)*scaling_factor
       end if
    else
    call write_log('Variable btractx_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'btracty',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%btracty)) then
       call write_log('  Loading btracty')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%btracty(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling btracty",GM_DIAGNOSTIC)
          data%stress%btracty(:,:) = &
               data%stress%btracty(:,:)*scaling_factor
       end if
    else
    call write_log('Variable btracty was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'btracty_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%btracty_extend)) then
       call write_log('  Loading btracty_extend')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%btracty_extend(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling btracty_extend",GM_DIAGNOSTIC)
          data%stress%btracty_extend(:,:) = &
               data%stress%btracty_extend(:,:)*scaling_factor
       end if
    else
    call write_log('Variable btracty_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bwat',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_hydro%bwat)) then
       call write_log('  Loading bwat')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_hydro%bwat, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bwat",GM_DIAGNOSTIC)
          data%basal_hydro%bwat = &
               data%basal_hydro%bwat*scaling_factor
       end if
    else
    call write_log('Variable bwat was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'bwat_diag',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_hydro%bwat_diag)) then
       call write_log('  Loading bwat_diag')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_hydro%bwat_diag, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling bwat_diag",GM_DIAGNOSTIC)
          data%basal_hydro%bwat_diag = &
               data%basal_hydro%bwat_diag*scaling_factor
       end if
    else
    call write_log('Variable bwat_diag was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'c_space_factor',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_physics%c_space_factor)) then
       call write_log('  Loading c_space_factor')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_physics%c_space_factor, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling c_space_factor",GM_DIAGNOSTIC)
          data%basal_physics%c_space_factor = &
               data%basal_physics%c_space_factor*scaling_factor
       end if
    else
    call write_log('Variable c_space_factor was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'calving_mask',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%calving%calving_mask)) then
       call write_log('  Loading calving_mask')
       status = distributed_get_var(NCI%id, varid, &
            data%calving%calving_mask, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling calving_mask",GM_DIAGNOSTIC)
          data%calving%calving_mask = &
               data%calving%calving_mask*scaling_factor
       end if
    else
    call write_log('Variable calving_mask was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'cell_area',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%cell_area)) then
       call write_log('  Loading cell_area')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%cell_area, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(len0*len0)
       else
          scaling_factor = scaling_factor/(len0*len0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling cell_area",GM_DIAGNOSTIC)
          data%geometry%cell_area = &
               data%geometry%cell_area*scaling_factor
       end if
    else
    call write_log('Variable cell_area was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'coulomb_c',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_physics%coulomb_c)) then
       call write_log('  Loading coulomb_c')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_physics%coulomb_c, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling coulomb_c",GM_DIAGNOSTIC)
          data%basal_physics%coulomb_c = &
               data%basal_physics%coulomb_c*scaling_factor
       end if
    else
    call write_log('Variable coulomb_c was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'coulomb_c_relax',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_physics%coulomb_c_relax)) then
       call write_log('  Loading coulomb_c_relax')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_physics%coulomb_c_relax, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling coulomb_c_relax",GM_DIAGNOSTIC)
          data%basal_physics%coulomb_c_relax = &
               data%basal_physics%coulomb_c_relax*scaling_factor
       end if
    else
    call write_log('Variable coulomb_c_relax was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'damage',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%calving%damage)) then
       call write_log('  Loading damage')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%calving%damage(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling damage",GM_DIAGNOSTIC)
             data%calving%damage(up,:,:) = &
                  data%calving%damage(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable damage was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'deltaT_ocn',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%ocean_data%deltaT_ocn)) then
       call write_log('  Loading deltaT_ocn')
       status = distributed_get_var(NCI%id, varid, &
            data%ocean_data%deltaT_ocn, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling deltaT_ocn",GM_DIAGNOSTIC)
          data%ocean_data%deltaT_ocn = &
               data%ocean_data%deltaT_ocn*scaling_factor
       end if
    else
    call write_log('Variable deltaT_ocn was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'dissip',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%dissip)) then
       call write_log('  Loading dissip')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%dissip(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scyr)
          else
             scaling_factor = scaling_factor/(scyr)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling dissip",GM_DIAGNOSTIC)
             data%temper%dissip(up,:,:) = &
                  data%temper%dissip(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable dissip was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'dissipstag',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%dissip)) then
       call write_log('  Loading dissipstag')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%dissip(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scyr)
          else
             scaling_factor = scaling_factor/(scyr)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling dissipstag",GM_DIAGNOSTIC)
             data%temper%dissip(up,:,:) = &
                  data%temper%dissip(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable dissipstag was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'divu',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%divu)) then
       call write_log('  Loading divu')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%divu, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scyr)
       else
          scaling_factor = scaling_factor/(scyr)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling divu",GM_DIAGNOSTIC)
          data%velocity%divu = &
               data%velocity%divu*scaling_factor
       end if
    else
    call write_log('Variable divu was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'dthck_dt',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%dthck_dt)) then
       call write_log('  Loading dthck_dt')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%dthck_dt, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scyr)
       else
          scaling_factor = scaling_factor/(scyr)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling dthck_dt",GM_DIAGNOSTIC)
          data%geometry%dthck_dt = &
               data%geometry%dthck_dt*scaling_factor
       end if
    else
    call write_log('Variable dthck_dt was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'dthck_dt_obs',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%dthck_dt_obs)) then
       call write_log('  Loading dthck_dt_obs')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%dthck_dt_obs, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling dthck_dt_obs",GM_DIAGNOSTIC)
          data%geometry%dthck_dt_obs = &
               data%geometry%dthck_dt_obs*scaling_factor
       end if
    else
    call write_log('Variable dthck_dt_obs was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'dthck_dt_obs_basin',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%dthck_dt_obs_basin)) then
       call write_log('  Loading dthck_dt_obs_basin')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%dthck_dt_obs_basin, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling dthck_dt_obs_basin",GM_DIAGNOSTIC)
          data%geometry%dthck_dt_obs_basin = &
               data%geometry%dthck_dt_obs_basin*scaling_factor
       end if
    else
    call write_log('Variable dthck_dt_obs_basin was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'efvs',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%efvs)) then
       call write_log('  Loading efvs')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%stress%efvs(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_efvs)
          else
             scaling_factor = scaling_factor/(scale_efvs)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling efvs",GM_DIAGNOSTIC)
             data%stress%efvs(up,:,:) = &
                  data%stress%efvs(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable efvs was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'f_effecpress_bwat',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_physics%f_effecpress_bwat)) then
       call write_log('  Loading f_effecpress_bwat')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_physics%f_effecpress_bwat, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling f_effecpress_bwat",GM_DIAGNOSTIC)
          data%basal_physics%f_effecpress_bwat = &
               data%basal_physics%f_effecpress_bwat*scaling_factor
       end if
    else
    call write_log('Variable f_effecpress_bwat was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'f_effecpress_bwat_target',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_physics%f_effecpress_bwat_target)) then
       call write_log('  Loading f_effecpress_bwat_target')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_physics%f_effecpress_bwat_target, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling f_effecpress_bwat_target",GM_DIAGNOSTIC)
          data%basal_physics%f_effecpress_bwat_target = &
               data%basal_physics%f_effecpress_bwat_target*scaling_factor
       end if
    else
    call write_log('Variable f_effecpress_bwat_target was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'f_effecpress_bwatflx',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_physics%f_effecpress_bwatflx)) then
       call write_log('  Loading f_effecpress_bwatflx')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_physics%f_effecpress_bwatflx, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling f_effecpress_bwatflx",GM_DIAGNOSTIC)
          data%basal_physics%f_effecpress_bwatflx = &
               data%basal_physics%f_effecpress_bwatflx*scaling_factor
       end if
    else
    call write_log('Variable f_effecpress_bwatflx was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'f_effecpress_ocean_p',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_physics%f_effecpress_ocean_p)) then
       call write_log('  Loading f_effecpress_ocean_p')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_physics%f_effecpress_ocean_p, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling f_effecpress_ocean_p",GM_DIAGNOSTIC)
          data%basal_physics%f_effecpress_ocean_p = &
               data%basal_physics%f_effecpress_ocean_p*scaling_factor
       end if
    else
    call write_log('Variable f_effecpress_ocean_p was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'ff_invert_mask',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%inversion%ff_invert_mask)) then
       call write_log('  Loading ff_invert_mask')
       status = distributed_get_var(NCI%id, varid, &
            data%inversion%ff_invert_mask, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling ff_invert_mask",GM_DIAGNOSTIC)
          data%inversion%ff_invert_mask = &
               data%inversion%ff_invert_mask*scaling_factor
       end if
    else
    call write_log('Variable ff_invert_mask was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'floating_thck_target',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%inversion%floating_thck_target)) then
       call write_log('  Loading floating_thck_target')
       status = distributed_get_var(NCI%id, varid, &
            data%inversion%floating_thck_target, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling floating_thck_target",GM_DIAGNOSTIC)
          data%inversion%floating_thck_target = &
               data%inversion%floating_thck_target*scaling_factor
       end if
    else
    call write_log('Variable floating_thck_target was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'flow_enhancement_factor',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%flow_enhancement_factor)) then
       call write_log('  Loading flow_enhancement_factor')
       status = distributed_get_var(NCI%id, varid, &
            data%temper%flow_enhancement_factor, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling flow_enhancement_factor",GM_DIAGNOSTIC)
          data%temper%flow_enhancement_factor = &
               data%temper%flow_enhancement_factor*scaling_factor
       end if
    else
    call write_log('Variable flow_enhancement_factor was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'flwa',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%flwa)) then
       call write_log('  Loading flwa')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%flwa(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_flwa)
          else
             scaling_factor = scaling_factor/(scale_flwa)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling flwa",GM_DIAGNOSTIC)
             data%temper%flwa(up,:,:) = &
                  data%temper%flwa(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable flwa was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'flwastag',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%flwa)) then
       call write_log('  Loading flwastag')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%flwa(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_flwa)
          else
             scaling_factor = scaling_factor/(scale_flwa)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling flwastag",GM_DIAGNOSTIC)
             data%temper%flwa(up,:,:) = &
                  data%temper%flwa(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable flwastag was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'ice_age',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%ice_age)) then
       call write_log('  Loading ice_age')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%geometry%ice_age(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(tim0/scyr)
          else
             scaling_factor = scaling_factor/(tim0/scyr)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling ice_age",GM_DIAGNOSTIC)
             data%geometry%ice_age(up,:,:) = &
                  data%geometry%ice_age(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable ice_age was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'ice_domain_mask',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%general%ice_domain_mask)) then
       call write_log('  Loading ice_domain_mask')
       status = distributed_get_var(NCI%id, varid, &
            data%general%ice_domain_mask, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling ice_domain_mask",GM_DIAGNOSTIC)
          data%general%ice_domain_mask = &
               data%general%ice_domain_mask*scaling_factor
       end if
    else
    call write_log('Variable ice_domain_mask was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'ice_fraction_retreat_mask',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%ice_fraction_retreat_mask)) then
       call write_log('  Loading ice_fraction_retreat_mask')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%ice_fraction_retreat_mask, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling ice_fraction_retreat_mask",GM_DIAGNOSTIC)
          data%geometry%ice_fraction_retreat_mask = &
               data%geometry%ice_fraction_retreat_mask*scaling_factor
       end if
    else
    call write_log('Variable ice_fraction_retreat_mask was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'kinbcmask',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%kinbcmask)) then
       call write_log('  Loading kinbcmask')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%kinbcmask(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling kinbcmask",GM_DIAGNOSTIC)
          data%velocity%kinbcmask(:,:) = &
               data%velocity%kinbcmask(:,:)*scaling_factor
       end if
    else
    call write_log('Variable kinbcmask was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'lat',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%general%lat)) then
       call write_log('  Loading lat')
       status = distributed_get_var(NCI%id, varid, &
            data%general%lat, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling lat",GM_DIAGNOSTIC)
          data%general%lat = &
               data%general%lat*scaling_factor
       end if
    else
    call write_log('Variable lat was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'litho_temp',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%lithot%temp)) then
       call write_log('  Loading litho_temp')
       status = distributed_get_var(NCI%id, varid, &
            data%lithot%temp, parallel, (/1,1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling litho_temp",GM_DIAGNOSTIC)
          data%lithot%temp = &
               data%lithot%temp*scaling_factor
       end if
    else
    call write_log('Variable litho_temp was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'load',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%isostasy%load)) then
       call write_log('  Loading load')
       status = distributed_get_var(NCI%id, varid, &
            data%isostasy%load, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling load",GM_DIAGNOSTIC)
          data%isostasy%load = &
               data%isostasy%load*scaling_factor
       end if
    else
    call write_log('Variable load was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'lon',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%general%lon)) then
       call write_log('  Loading lon')
       status = distributed_get_var(NCI%id, varid, &
            data%general%lon, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling lon",GM_DIAGNOSTIC)
          data%general%lon = &
               data%general%lon*scaling_factor
       end if
    else
    call write_log('Variable lon was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'marine_connection_mask_isolated',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%marine_connection_mask_isolated)) then
       call write_log('  Loading marine_connection_mask_isolated')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%marine_connection_mask_isolated, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling marine_connection_mask_isolated",GM_DIAGNOSTIC)
          data%geometry%marine_connection_mask_isolated = &
               data%geometry%marine_connection_mask_isolated*scaling_factor
       end if
    else
    call write_log('Variable marine_connection_mask_isolated was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'overwrite_acab_mask',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%overwrite_acab_mask)) then
       call write_log('  Loading overwrite_acab_mask')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%overwrite_acab_mask, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling overwrite_acab_mask",GM_DIAGNOSTIC)
          data%climate%overwrite_acab_mask = &
               data%climate%overwrite_acab_mask*scaling_factor
       end if
    else
    call write_log('Variable overwrite_acab_mask was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'powerlaw_c',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_physics%powerlaw_c)) then
       call write_log('  Loading powerlaw_c')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_physics%powerlaw_c, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling powerlaw_c",GM_DIAGNOSTIC)
          data%basal_physics%powerlaw_c = &
               data%basal_physics%powerlaw_c*scaling_factor
       end if
    else
    call write_log('Variable powerlaw_c was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'reference_thck',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%reference_thck)) then
       call write_log('  Loading reference_thck')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%reference_thck, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling reference_thck",GM_DIAGNOSTIC)
          data%geometry%reference_thck = &
               data%geometry%reference_thck*scaling_factor
       end if
    else
    call write_log('Variable reference_thck was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'relx',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%isostasy%relx)) then
       call write_log('  Loading relx')
       status = distributed_get_var(NCI%id, varid, &
            data%isostasy%relx, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling relx",GM_DIAGNOSTIC)
          data%isostasy%relx = &
               data%isostasy%relx*scaling_factor
       end if
    else
    call write_log('Variable relx was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'smb',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%smb)) then
       call write_log('  Loading smb')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%smb, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(1.0)
       else
          scaling_factor = scaling_factor/(1.0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling smb",GM_DIAGNOSTIC)
          data%climate%smb = &
               data%climate%smb*scaling_factor
       end if
    else
    call write_log('Variable smb was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'smb_3d',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%smb_3d)) then
       call write_log('  Loading smb_3d')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%smb_3d(up,:,:), parallel, (/1,1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(1.0)
       else
          scaling_factor = scaling_factor/(1.0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling smb_3d",GM_DIAGNOSTIC)
          data%climate%smb_3d(up,:,:) = &
               data%climate%smb_3d(up,:,:)*scaling_factor
       end if
    else
    call write_log('Variable smb_3d was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'smb_anomaly',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%smb_anomaly)) then
       call write_log('  Loading smb_anomaly')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%smb_anomaly, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(1.0)
       else
          scaling_factor = scaling_factor/(1.0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling smb_anomaly",GM_DIAGNOSTIC)
          data%climate%smb_anomaly = &
               data%climate%smb_anomaly*scaling_factor
       end if
    else
    call write_log('Variable smb_anomaly was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'smb_gradz',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%smb_gradz)) then
       call write_log('  Loading smb_gradz')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%smb_gradz, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(1.0/thk0)
       else
          scaling_factor = scaling_factor/(1.0/thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling smb_gradz",GM_DIAGNOSTIC)
          data%climate%smb_gradz = &
               data%climate%smb_gradz*scaling_factor
       end if
    else
    call write_log('Variable smb_gradz was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'smb_levels',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%smb_levels)) then
       call write_log('  Loading smb_levels')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%smb_levels, parallel, (/1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling smb_levels",GM_DIAGNOSTIC)
          data%climate%smb_levels = &
               data%climate%smb_levels*scaling_factor
       end if
    else
    call write_log('Variable smb_levels was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'smb_ref',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%smb_ref)) then
       call write_log('  Loading smb_ref')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%smb_ref, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(1.0)
       else
          scaling_factor = scaling_factor/(1.0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling smb_ref",GM_DIAGNOSTIC)
          data%climate%smb_ref = &
               data%climate%smb_ref*scaling_factor
       end if
    else
    call write_log('Variable smb_ref was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'smb_reference_usrf',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%climate%smb_reference_usrf)) then
       call write_log('  Loading smb_reference_usrf')
       status = distributed_get_var(NCI%id, varid, &
            data%climate%smb_reference_usrf, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling smb_reference_usrf",GM_DIAGNOSTIC)
          data%climate%smb_reference_usrf = &
               data%climate%smb_reference_usrf*scaling_factor
       end if
    else
    call write_log('Variable smb_reference_usrf was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'soft',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%bed_softness)) then
       call write_log('  Loading soft')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%bed_softness, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_btrc)
       else
          scaling_factor = scaling_factor/(scale_btrc)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling soft",GM_DIAGNOSTIC)
          data%velocity%bed_softness = &
               data%velocity%bed_softness*scaling_factor
       end if
    else
    call write_log('Variable soft was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'tau_eigen1',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%calving%tau_eigen1)) then
       call write_log('  Loading tau_eigen1')
       status = distributed_get_var(NCI%id, varid, &
            data%calving%tau_eigen1, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling tau_eigen1",GM_DIAGNOSTIC)
          data%calving%tau_eigen1 = &
               data%calving%tau_eigen1*scaling_factor
       end if
    else
    call write_log('Variable tau_eigen1 was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'tau_eigen2',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%calving%tau_eigen2)) then
       call write_log('  Loading tau_eigen2')
       status = distributed_get_var(NCI%id, varid, &
            data%calving%tau_eigen2, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling tau_eigen2",GM_DIAGNOSTIC)
          data%calving%tau_eigen2 = &
               data%calving%tau_eigen2*scaling_factor
       end if
    else
    call write_log('Variable tau_eigen2 was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'taudx',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%taudx)) then
       call write_log('  Loading taudx')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%taudx(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling taudx",GM_DIAGNOSTIC)
          data%stress%taudx(:,:) = &
               data%stress%taudx(:,:)*scaling_factor
       end if
    else
    call write_log('Variable taudx was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'taudy',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%stress%taudy)) then
       call write_log('  Loading taudy')
       status = distributed_get_var(NCI%id, varid, &
            data%stress%taudy(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling taudy",GM_DIAGNOSTIC)
          data%stress%taudy(:,:) = &
               data%stress%taudy(:,:)*scaling_factor
       end if
    else
    call write_log('Variable taudy was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'tauf',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_physics%mintauf)) then
       call write_log('  Loading tauf')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_physics%mintauf, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_tau)
       else
          scaling_factor = scaling_factor/(scale_tau)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling tauf",GM_DIAGNOSTIC)
          data%basal_physics%mintauf = &
               data%basal_physics%mintauf*scaling_factor
       end if
    else
    call write_log('Variable tauf was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'temp',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%temp)) then
       call write_log('  Loading temp')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%temp(up,1:data%general%ewn,1:data%general%nsn), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling temp",GM_DIAGNOSTIC)
             data%temper%temp(up,1:data%general%ewn,1:data%general%nsn) = &
                  data%temper%temp(up,1:data%general%ewn,1:data%general%nsn)*scaling_factor
          end if
       end do
    else
    call write_log('Variable temp was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'tempstag',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%temp)) then
       call write_log('  Loading tempstag')
       do up=0,NCI%nstagwbndlevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%temp(up,1:data%general%ewn,1:data%general%nsn), parallel, (/1,1,up+1,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling tempstag",GM_DIAGNOSTIC)
             data%temper%temp(up,1:data%general%ewn,1:data%general%nsn) = &
                  data%temper%temp(up,1:data%general%ewn,1:data%general%nsn)*scaling_factor
          end if
       end do
    else
    call write_log('Variable tempstag was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'tempunstag',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%tempunstag)) then
       call write_log('  Loading tempunstag')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%tempunstag(up,1:data%general%ewn,1:data%general%nsn), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling tempunstag",GM_DIAGNOSTIC)
             data%temper%tempunstag(up,1:data%general%ewn,1:data%general%nsn) = &
                  data%temper%tempunstag(up,1:data%general%ewn,1:data%general%nsn)*scaling_factor
          end if
       end do
    else
    call write_log('Variable tempunstag was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'thck_calving_threshold',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%calving%thck_calving_threshold)) then
       call write_log('  Loading thck_calving_threshold')
       status = distributed_get_var(NCI%id, varid, &
            data%calving%thck_calving_threshold, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling thck_calving_threshold",GM_DIAGNOSTIC)
          data%calving%thck_calving_threshold = &
               data%calving%thck_calving_threshold*scaling_factor
       end if
    else
    call write_log('Variable thck_calving_threshold was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'thermal_forcing',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%ocean_data%thermal_forcing)) then
       call write_log('  Loading thermal_forcing')
       do up=1,NCI%nzocn
          status = distributed_get_var(NCI%id, varid, &
               data%ocean_data%thermal_forcing(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling thermal_forcing",GM_DIAGNOSTIC)
             data%ocean_data%thermal_forcing(up,:,:) = &
                  data%ocean_data%thermal_forcing(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable thermal_forcing was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'thk',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%thck)) then
       call write_log('  Loading thk')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%thck, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling thk",GM_DIAGNOSTIC)
          data%geometry%thck = &
               data%geometry%thck*scaling_factor
       end if
    else
    call write_log('Variable thk was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'thkmask',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%thkmask)) then
       call write_log('  Loading thkmask')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%thkmask, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling thkmask",GM_DIAGNOSTIC)
          data%geometry%thkmask = &
               data%geometry%thkmask*scaling_factor
       end if
    else
    call write_log('Variable thkmask was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'topg',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%topg)) then
       call write_log('  Loading topg')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%topg, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling topg",GM_DIAGNOSTIC)
          data%geometry%topg = &
               data%geometry%topg*scaling_factor
       end if
    else
    call write_log('Variable topg was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'topg_stdev',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%topg_stdev)) then
       call write_log('  Loading topg_stdev')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%topg_stdev, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling topg_stdev",GM_DIAGNOSTIC)
          data%geometry%topg_stdev = &
               data%geometry%topg_stdev*scaling_factor
       end if
    else
    call write_log('Variable topg_stdev was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'unstagbeta',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%unstagbeta)) then
       call write_log('  Loading unstagbeta')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%unstagbeta, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_beta)
       else
          scaling_factor = scaling_factor/(scale_beta)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling unstagbeta",GM_DIAGNOSTIC)
          data%velocity%unstagbeta = &
               data%velocity%unstagbeta*scaling_factor
       end if
    else
    call write_log('Variable unstagbeta was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'usfc_obs',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%usfc_obs)) then
       call write_log('  Loading usfc_obs')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%usfc_obs, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling usfc_obs",GM_DIAGNOSTIC)
          data%velocity%usfc_obs = &
               data%velocity%usfc_obs*scaling_factor
       end if
    else
    call write_log('Variable usfc_obs was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'usrf_obs',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%usrf_obs)) then
       call write_log('  Loading usrf_obs')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%usrf_obs, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling usrf_obs",GM_DIAGNOSTIC)
          data%geometry%usrf_obs = &
               data%geometry%usrf_obs*scaling_factor
       end if
    else
    call write_log('Variable usrf_obs was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'usurf',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%geometry%usrf)) then
       call write_log('  Loading usurf')
       status = distributed_get_var(NCI%id, varid, &
            data%geometry%usrf, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(thk0)
       else
          scaling_factor = scaling_factor/(thk0)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling usurf",GM_DIAGNOSTIC)
          data%geometry%usrf = &
               data%geometry%usrf*scaling_factor
       end if
    else
    call write_log('Variable usurf was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'uvel',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%uvel)) then
       call write_log('  Loading uvel')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%uvel(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_uvel)
          else
             scaling_factor = scaling_factor/(scale_uvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling uvel",GM_DIAGNOSTIC)
             data%velocity%uvel(up,:,:) = &
                  data%velocity%uvel(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable uvel was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'uvel_2d',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%uvel_2d)) then
       call write_log('  Loading uvel_2d')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%uvel_2d(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling uvel_2d",GM_DIAGNOSTIC)
          data%velocity%uvel_2d(:,:) = &
               data%velocity%uvel_2d(:,:)*scaling_factor
       end if
    else
    call write_log('Variable uvel_2d was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'uvel_2d_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%uvel_2d_extend)) then
       call write_log('  Loading uvel_2d_extend')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%uvel_2d_extend(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling uvel_2d_extend",GM_DIAGNOSTIC)
          data%velocity%uvel_2d_extend(:,:) = &
               data%velocity%uvel_2d_extend(:,:)*scaling_factor
       end if
    else
    call write_log('Variable uvel_2d_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'uvel_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%uvel_extend)) then
       call write_log('  Loading uvel_extend')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%uvel_extend(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_uvel)
          else
             scaling_factor = scaling_factor/(scale_uvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling uvel_extend",GM_DIAGNOSTIC)
             data%velocity%uvel_extend(up,:,:) = &
                  data%velocity%uvel_extend(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable uvel_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'velo_sfc_obs',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%velo_sfc_obs)) then
       call write_log('  Loading velo_sfc_obs')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%velo_sfc_obs, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling velo_sfc_obs",GM_DIAGNOSTIC)
          data%velocity%velo_sfc_obs = &
               data%velocity%velo_sfc_obs*scaling_factor
       end if
    else
    call write_log('Variable velo_sfc_obs was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'vsfc_obs',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%vsfc_obs)) then
       call write_log('  Loading vsfc_obs')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%vsfc_obs, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling vsfc_obs",GM_DIAGNOSTIC)
          data%velocity%vsfc_obs = &
               data%velocity%vsfc_obs*scaling_factor
       end if
    else
    call write_log('Variable vsfc_obs was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'vvel',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%vvel)) then
       call write_log('  Loading vvel')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%vvel(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_uvel)
          else
             scaling_factor = scaling_factor/(scale_uvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling vvel",GM_DIAGNOSTIC)
             data%velocity%vvel(up,:,:) = &
                  data%velocity%vvel(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable vvel was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'vvel_2d',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%vvel_2d)) then
       call write_log('  Loading vvel_2d')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%vvel_2d(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling vvel_2d",GM_DIAGNOSTIC)
          data%velocity%vvel_2d(:,:) = &
               data%velocity%vvel_2d(:,:)*scaling_factor
       end if
    else
    call write_log('Variable vvel_2d was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'vvel_2d_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%vvel_2d_extend)) then
       call write_log('  Loading vvel_2d_extend')
       status = distributed_get_var(NCI%id, varid, &
            data%velocity%vvel_2d_extend(:,:), parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0/(scale_uvel)
       else
          scaling_factor = scaling_factor/(scale_uvel)
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling vvel_2d_extend",GM_DIAGNOSTIC)
          data%velocity%vvel_2d_extend(:,:) = &
               data%velocity%vvel_2d_extend(:,:)*scaling_factor
       end if
    else
    call write_log('Variable vvel_2d_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'vvel_extend',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%vvel_extend)) then
       call write_log('  Loading vvel_extend')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%vvel_extend(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_uvel)
          else
             scaling_factor = scaling_factor/(scale_uvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling vvel_extend",GM_DIAGNOSTIC)
             data%velocity%vvel_extend(up,:,:) = &
                  data%velocity%vvel_extend(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable vvel_extend was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'warm_ocean_mask',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%basal_melt%warm_ocean_mask)) then
       call write_log('  Loading warm_ocean_mask')
       status = distributed_get_var(NCI%id, varid, &
            data%basal_melt%warm_ocean_mask, parallel, (/1,1,infile%current_time/))
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
       if (status.ne.NF90_NOERR) then
          scaling_factor = 1.0d0
       end if
       if (abs(scaling_factor-1.0d0).gt.1.d-17) then
          call write_log("scaling warm_ocean_mask",GM_DIAGNOSTIC)
          data%basal_melt%warm_ocean_mask = &
               data%basal_melt%warm_ocean_mask*scaling_factor
       end if
    else
    call write_log('Variable warm_ocean_mask was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'waterfrac',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%temper%waterfrac)) then
       call write_log('  Loading waterfrac')
       do up=1,NCI%nstaglevel
          status = distributed_get_var(NCI%id, varid, &
               data%temper%waterfrac(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling waterfrac",GM_DIAGNOSTIC)
             data%temper%waterfrac(up,:,:) = &
                  data%temper%waterfrac(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable waterfrac was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'wgrd',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%wgrd)) then
       call write_log('  Loading wgrd')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%wgrd(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_wvel)
          else
             scaling_factor = scaling_factor/(scale_wvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling wgrd",GM_DIAGNOSTIC)
             data%velocity%wgrd(up,:,:) = &
                  data%velocity%wgrd(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable wgrd was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if

    status = parallel_inq_varid(NCI%id,'wvel',varid)
    if (status .eq. nf90_noerr) then
    if (is_enabled(data%velocity%wvel)) then
       call write_log('  Loading wvel')
       do up=1,NCI%nlevel
          status = distributed_get_var(NCI%id, varid, &
               data%velocity%wvel(up,:,:), parallel, (/1,1,up,infile%current_time/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_get_att(NCI%id, varid,'scale_factor',scaling_factor)
          if (status.ne.NF90_NOERR) then
             scaling_factor = 1.0d0/(scale_wvel)
          else
             scaling_factor = scaling_factor/(scale_wvel)
          end if
          if (abs(scaling_factor-1.0d0).gt.1.d-17) then
             call write_log("scaling wvel",GM_DIAGNOSTIC)
             data%velocity%wvel(up,:,:) = &
                  data%velocity%wvel(up,:,:)*scaling_factor
          end if
       end do
    else
    call write_log('Variable wvel was specified for input but it is &
         &inappropriate for your config settings.  It will be excluded from the input.', GM_WARNING)
    end if

    end if


  end subroutine glide_io_read

  subroutine glide_io_checkdim(infile,model,data)

    ! check if dimension sizes in file match dims of model
    use cism_parallel, only: parallel_type, parallel_inq_dimid, parallel_inquire_dimension
    use glimmer_log
    implicit none
    type(glimmer_nc_input), pointer :: infile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: model
    type(glide_global_type), optional :: data

    integer status,dimid,dimsize
    character(len=150) message

    ! check dimensions
    status = parallel_inq_dimid(NCI%id,'level',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%general%upn) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size level does not match: ', &
               model%general%upn
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'lithoz',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%lithot%nlayer) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size lithoz does not match: ', &
               model%lithot%nlayer
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'nlev_smb',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%climate%nlev_smb) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size nlev_smb does not match: ', &
               model%climate%nlev_smb
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'staglevel',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%general%upn-1) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size staglevel does not match: ', &
               model%general%upn-1
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'stagwbndlevel',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%general%upn+1) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size stagwbndlevel does not match: ', &
               model%general%upn+1
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'x0',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%parallel%global_ewn-1) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size x0 does not match: ', &
               model%parallel%global_ewn-1
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'x1',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%parallel%global_ewn) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size x1 does not match: ', &
               model%parallel%global_ewn
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'y0',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%parallel%global_nsn-1) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size y0 does not match: ', &
               model%parallel%global_nsn-1
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'y1',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.model%parallel%global_nsn) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size y1 does not match: ', &
               model%parallel%global_nsn
          call write_log(message,GM_FATAL)
       end if
    end if
    status = parallel_inq_dimid(NCI%id,'zocn',dimid)
    if (dimid.gt.0) then
       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
       if (dimsize.ne.data%ocean_data%nzocn) then
          write(message,*) 'Error, reading file ',trim(NCI%filename),' size zocn does not match: ', &
               data%ocean_data%nzocn
          call write_log(message,GM_FATAL)
       end if
    end if
  end subroutine glide_io_checkdim

  !*****************************************************************************
  ! calculating time averages
  !*****************************************************************************  
#ifdef HAVE_AVG
  subroutine glide_avg_accumulate(outfile,data,model)

    ! Accumulate time averages for 'tavg' variables
    ! NOTE: This subroutine works for tavg variables that are written to exactly one output file.
    !       If we try to write tavg variables to multiple output files, the accumulated values
    !        will be too large, because this subroutine will be called more than once per time step.
    ! TODO: Write code to check for doubly listed tavg variables and throw a fatal error.

    use cism_parallel, only: parallel_inq_varid
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: model
    type(glide_global_type) :: data

    ! local variables
    real(dp) :: factor
    integer status, varid

    ! increase total time
    outfile%total_time = outfile%total_time + model%numerics%tinc
    factor = model%numerics%tinc

    ! accumulate acab_applied
    status = parallel_inq_varid(NCO%id,'acab_applied_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%climate%acab_applied_tavg = data%climate%acab_applied_tavg + factor * data%climate%acab_applied
    end if

    ! accumulate basal_mbal_flux
    status = parallel_inq_varid(NCO%id,'basal_mbal_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%basal_mbal_flux_tavg = data%geometry%basal_mbal_flux_tavg + factor * data%geometry%basal_mbal_flux
    end if

    ! accumulate bmlt_applied
    status = parallel_inq_varid(NCO%id,'bmlt_applied_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%basal_melt%bmlt_applied_tavg = data%basal_melt%bmlt_applied_tavg + factor * data%basal_melt%bmlt_applied
    end if

    ! accumulate calving_flux
    status = parallel_inq_varid(NCO%id,'calving_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%calving_flux_tavg = data%geometry%calving_flux_tavg + factor * data%geometry%calving_flux
    end if

    ! accumulate calving_rate
    status = parallel_inq_varid(NCO%id,'calving_rate_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%calving%calving_rate_tavg = data%calving%calving_rate_tavg + factor * data%calving%calving_rate
    end if

    ! accumulate gl_flux
    status = parallel_inq_varid(NCO%id,'gl_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%gl_flux_tavg = data%geometry%gl_flux_tavg + factor * data%geometry%gl_flux
    end if

    ! accumulate sfc_mbal_flux
    status = parallel_inq_varid(NCO%id,'sfc_mbal_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%sfc_mbal_flux_tavg = data%geometry%sfc_mbal_flux_tavg + factor * data%geometry%sfc_mbal_flux
    end if

    ! accumulate total_bmb_flux
    status = parallel_inq_varid(NCO%id,'total_bmb_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%total_bmb_flux_tavg = data%geometry%total_bmb_flux_tavg + factor * data%geometry%total_bmb_flux
    end if

    ! accumulate total_calving_flux
    status = parallel_inq_varid(NCO%id,'total_calving_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%total_calving_flux_tavg = data%geometry%total_calving_flux_tavg + factor * data%geometry%total_calving_flux
    end if

    ! accumulate total_gl_flux
    status = parallel_inq_varid(NCO%id,'total_gl_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%total_gl_flux_tavg = data%geometry%total_gl_flux_tavg + factor * data%geometry%total_gl_flux
    end if

    ! accumulate total_smb_flux
    status = parallel_inq_varid(NCO%id,'total_smb_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%total_smb_flux_tavg = data%geometry%total_smb_flux_tavg + factor * data%geometry%total_smb_flux
    end if

  end subroutine glide_avg_accumulate

  subroutine glide_avg_reset(outfile,data)

    use cism_parallel, only: parallel_inq_varid
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    ! structure containg output netCDF descriptor
    type(glide_global_type) :: data

    ! local variables
    integer status, varid

    ! reset total time
    outfile%total_time = 0.d0

    ! reset acab_applied
    status = parallel_inq_varid(NCO%id,'acab_applied_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%climate%acab_applied_tavg = 0.
    end if

    ! reset basal_mbal_flux
    status = parallel_inq_varid(NCO%id,'basal_mbal_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%basal_mbal_flux_tavg = 0.
    end if

    ! reset bmlt_applied
    status = parallel_inq_varid(NCO%id,'bmlt_applied_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%basal_melt%bmlt_applied_tavg = 0.
    end if

    ! reset calving_flux
    status = parallel_inq_varid(NCO%id,'calving_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%calving_flux_tavg = 0.
    end if

    ! reset calving_rate
    status = parallel_inq_varid(NCO%id,'calving_rate_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%calving%calving_rate_tavg = 0.
    end if

    ! reset gl_flux
    status = parallel_inq_varid(NCO%id,'gl_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%gl_flux_tavg = 0.
    end if

    ! reset sfc_mbal_flux
    status = parallel_inq_varid(NCO%id,'sfc_mbal_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%sfc_mbal_flux_tavg = 0.
    end if

    ! reset total_bmb_flux
    status = parallel_inq_varid(NCO%id,'total_bmb_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%total_bmb_flux_tavg = 0.
    end if

    ! reset total_calving_flux
    status = parallel_inq_varid(NCO%id,'total_calving_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%total_calving_flux_tavg = 0.
    end if

    ! reset total_gl_flux
    status = parallel_inq_varid(NCO%id,'total_gl_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%total_gl_flux_tavg = 0.
    end if

    ! reset total_smb_flux
    status = parallel_inq_varid(NCO%id,'total_smb_flux_tavg',varid)
    if (status .eq. nf90_noerr) then
       data%geometry%total_smb_flux_tavg = 0.
    end if

  end subroutine glide_avg_reset
#endif

  !*********************************************************************
  ! some private procedures
  !*********************************************************************

  !> apply default type to be used in netCDF file
  integer function get_xtype(outfile,xtype)
    implicit none
    type(glimmer_nc_output), pointer :: outfile !< derived type holding information about output file
    integer, intent(in) :: xtype                !< the external netCDF type

    get_xtype = xtype
    
    if (xtype.eq.NF90_REAL .and. outfile%default_xtype.eq.NF90_DOUBLE) then
       get_xtype = NF90_DOUBLE
    end if
    if (xtype.eq.NF90_DOUBLE .and. outfile%default_xtype.eq.NF90_REAL) then
       get_xtype = NF90_REAL
    end if
  end function get_xtype

  !*********************************************************************
  ! lots of accessor subroutines follow
  !*********************************************************************
  subroutine glide_get_S_ambient(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%plume%S_ambient(1:data%general%ewn,1:data%general%nsn)
  end subroutine glide_get_S_ambient

  subroutine glide_get_T_ambient(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%plume%T_ambient(1:data%general%ewn,1:data%general%nsn)
  end subroutine glide_get_T_ambient

  subroutine glide_get_acab(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%climate%acab)
  end subroutine glide_get_acab

  subroutine glide_set_acab(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%acab = inarray/(scale_acab)
  end subroutine glide_set_acab

  subroutine glide_get_acab_anomaly(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%climate%acab_anomaly)
  end subroutine glide_get_acab_anomaly

  subroutine glide_set_acab_anomaly(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%acab_anomaly = inarray/(scale_acab)
  end subroutine glide_set_acab_anomaly

  subroutine glide_get_acab_applied(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%climate%acab_applied)
  end subroutine glide_get_acab_applied

  subroutine glide_set_acab_applied(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%acab_applied = inarray/(scale_acab)
  end subroutine glide_set_acab_applied

  subroutine glide_get_acab_corrected(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%climate%acab_corrected)
  end subroutine glide_get_acab_corrected

  subroutine glide_set_acab_corrected(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%acab_corrected = inarray/(scale_acab)
  end subroutine glide_set_acab_corrected

  subroutine glide_get_acab_gradz(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab/thk0)*(data%climate%acab_gradz)
  end subroutine glide_get_acab_gradz

  subroutine glide_set_acab_gradz(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%acab_gradz = inarray/(scale_acab/thk0)
  end subroutine glide_set_acab_gradz

  subroutine glide_get_acab_ref(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%climate%acab_ref)
  end subroutine glide_get_acab_ref

  subroutine glide_set_acab_ref(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%acab_ref = inarray/(scale_acab)
  end subroutine glide_set_acab_ref

  subroutine glide_get_adv_cfl_dt(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%numerics%adv_cfl_dt
  end subroutine glide_get_adv_cfl_dt

  subroutine glide_set_adv_cfl_dt(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%numerics%adv_cfl_dt = inarray
  end subroutine glide_set_adv_cfl_dt

  subroutine glide_get_area_factor(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%projection%stere%area_factor
  end subroutine glide_get_area_factor

  subroutine glide_set_area_factor(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%projection%stere%area_factor = inarray
  end subroutine glide_set_area_factor

  subroutine glide_get_artm(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%climate%artm
  end subroutine glide_get_artm

  subroutine glide_set_artm(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%artm = inarray
  end subroutine glide_set_artm

  subroutine glide_get_artm_anomaly(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%climate%artm_anomaly
  end subroutine glide_get_artm_anomaly

  subroutine glide_set_artm_anomaly(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%artm_anomaly = inarray
  end subroutine glide_set_artm_anomaly

  subroutine glide_get_artm_gradz(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1./thk0)*(data%climate%artm_gradz)
  end subroutine glide_get_artm_gradz

  subroutine glide_set_artm_gradz(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%artm_gradz = inarray/(1./thk0)
  end subroutine glide_set_artm_gradz

  subroutine glide_get_artm_ref(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%climate%artm_ref
  end subroutine glide_get_artm_ref

  subroutine glide_set_artm_ref(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%artm_ref = inarray
  end subroutine glide_set_artm_ref

  subroutine glide_get_basal_mbal_flux(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%basal_mbal_flux
  end subroutine glide_get_basal_mbal_flux

  subroutine glide_set_basal_mbal_flux(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%basal_mbal_flux = inarray
  end subroutine glide_set_basal_mbal_flux

  subroutine glide_get_basin_multiplier_array(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_melt%basin_multiplier_array
  end subroutine glide_get_basin_multiplier_array

  subroutine glide_set_basin_multiplier_array(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_melt%basin_multiplier_array = inarray
  end subroutine glide_set_basin_multiplier_array

  subroutine glide_get_basin_number(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%ocean_data%basin_number
  end subroutine glide_get_basin_number

  subroutine glide_get_beta(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_beta)*(data%velocity%beta)
  end subroutine glide_get_beta

  subroutine glide_set_beta(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%beta = inarray/(scale_beta)
  end subroutine glide_set_beta

  subroutine glide_get_beta_internal(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_beta)*(data%velocity%beta_internal)
  end subroutine glide_get_beta_internal

  subroutine glide_set_beta_internal(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%beta_internal = inarray/(scale_beta)
  end subroutine glide_set_beta_internal

  subroutine glide_get_bfricflx(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%temper%bfricflx
  end subroutine glide_get_bfricflx

  subroutine glide_set_bfricflx(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%bfricflx = inarray
  end subroutine glide_set_bfricflx

  subroutine glide_get_bheatflx(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%temper%bheatflx
  end subroutine glide_get_bheatflx

  subroutine glide_set_bheatflx(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%bheatflx = inarray
  end subroutine glide_set_bheatflx

  subroutine glide_get_bmlt(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%basal_melt%bmlt)
  end subroutine glide_get_bmlt

  subroutine glide_set_bmlt(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_melt%bmlt = inarray/(scale_acab)
  end subroutine glide_set_bmlt

  subroutine glide_get_bmlt_applied(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%basal_melt%bmlt_applied)
  end subroutine glide_get_bmlt_applied

  subroutine glide_set_bmlt_applied(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_melt%bmlt_applied = inarray/(scale_acab)
  end subroutine glide_set_bmlt_applied

  subroutine glide_get_bmlt_float(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%basal_melt%bmlt_float)
  end subroutine glide_get_bmlt_float

  subroutine glide_set_bmlt_float(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_melt%bmlt_float = inarray/(scale_acab)
  end subroutine glide_set_bmlt_float

  subroutine glide_get_bmlt_float_anomaly(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%basal_melt%bmlt_float_anomaly)
  end subroutine glide_get_bmlt_float_anomaly

  subroutine glide_set_bmlt_float_anomaly(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_melt%bmlt_float_anomaly = inarray/(scale_acab)
  end subroutine glide_set_bmlt_float_anomaly

  subroutine glide_get_bmlt_float_external(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%basal_melt%bmlt_float_external)
  end subroutine glide_get_bmlt_float_external

  subroutine glide_set_bmlt_float_external(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_melt%bmlt_float_external = inarray/(scale_acab)
  end subroutine glide_set_bmlt_float_external

  subroutine glide_get_bmlt_ground(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%basal_melt%bmlt_ground)
  end subroutine glide_get_bmlt_ground

  subroutine glide_set_bmlt_ground(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_melt%bmlt_ground = inarray/(scale_acab)
  end subroutine glide_set_bmlt_ground

  subroutine glide_get_bpmp(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%temper%bpmp
  end subroutine glide_get_bpmp

  subroutine glide_set_bpmp(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%bpmp = inarray
  end subroutine glide_set_bpmp

  subroutine glide_get_btemp(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%temper%btemp
  end subroutine glide_get_btemp

  subroutine glide_set_btemp(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%btemp = inarray
  end subroutine glide_set_btemp

  subroutine glide_get_btemp_float(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%temper%btemp_float
  end subroutine glide_get_btemp_float

  subroutine glide_set_btemp_float(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%btemp_float = inarray
  end subroutine glide_set_btemp_float

  subroutine glide_get_btemp_ground(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%temper%btemp_ground
  end subroutine glide_get_btemp_ground

  subroutine glide_set_btemp_ground(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%btemp_ground = inarray
  end subroutine glide_set_btemp_ground

  subroutine glide_get_btract(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%btract(:,:))
  end subroutine glide_get_btract

  subroutine glide_set_btract(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%btract(:,:) = inarray/(scale_tau)
  end subroutine glide_set_btract

  subroutine glide_get_btractx(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%btractx(:,:))
  end subroutine glide_get_btractx

  subroutine glide_set_btractx(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%btractx(:,:) = inarray/(scale_tau)
  end subroutine glide_set_btractx

  subroutine glide_get_btractx_extend(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%btractx_extend(:,:))
  end subroutine glide_get_btractx_extend

  subroutine glide_set_btractx_extend(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%btractx_extend(:,:) = inarray/(scale_tau)
  end subroutine glide_set_btractx_extend

  subroutine glide_get_btracty(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%btracty(:,:))
  end subroutine glide_get_btracty

  subroutine glide_set_btracty(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%btracty(:,:) = inarray/(scale_tau)
  end subroutine glide_set_btracty

  subroutine glide_get_btracty_extend(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%btracty_extend(:,:))
  end subroutine glide_get_btracty_extend

  subroutine glide_set_btracty_extend(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%btracty_extend(:,:) = inarray/(scale_tau)
  end subroutine glide_set_btracty_extend

  subroutine glide_get_btrc(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_btrc)*(data%velocity%btrc)
  end subroutine glide_get_btrc

  subroutine glide_set_btrc(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%btrc = inarray/(scale_btrc)
  end subroutine glide_set_btrc

  subroutine glide_get_bwat(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%basal_hydro%bwat)
  end subroutine glide_get_bwat

  subroutine glide_set_bwat(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_hydro%bwat = inarray/(thk0)
  end subroutine glide_set_bwat

  subroutine glide_get_bwat_diag(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%basal_hydro%bwat_diag)
  end subroutine glide_get_bwat_diag

  subroutine glide_set_bwat_diag(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_hydro%bwat_diag = inarray/(thk0)
  end subroutine glide_set_bwat_diag

  subroutine glide_get_bwatflx(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_hydro%bwatflx
  end subroutine glide_get_bwatflx

  subroutine glide_set_bwatflx(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_hydro%bwatflx = inarray
  end subroutine glide_set_bwatflx

  subroutine glide_get_c_flux_array(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_hydro%c_flux_array
  end subroutine glide_get_c_flux_array

  subroutine glide_set_c_flux_array(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_hydro%c_flux_array = inarray
  end subroutine glide_set_c_flux_array

  subroutine glide_get_c_space_factor(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_physics%c_space_factor
  end subroutine glide_get_c_space_factor

  subroutine glide_set_c_space_factor(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%c_space_factor = inarray
  end subroutine glide_set_c_space_factor

  subroutine glide_get_calving_flux(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%calving_flux
  end subroutine glide_get_calving_flux

  subroutine glide_set_calving_flux(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%calving_flux = inarray
  end subroutine glide_set_calving_flux

  subroutine glide_get_calving_lateral(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scyr)*(data%calving%lateral_rate)
  end subroutine glide_get_calving_lateral

  subroutine glide_set_calving_lateral(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%calving%lateral_rate = inarray/(scyr)
  end subroutine glide_set_calving_lateral

  subroutine glide_get_calving_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%calving%calving_mask
  end subroutine glide_get_calving_mask

  subroutine glide_set_calving_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%calving%calving_mask = inarray
  end subroutine glide_set_calving_mask

  subroutine glide_get_calving_rate(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%calving%calving_rate
  end subroutine glide_get_calving_rate

  subroutine glide_set_calving_rate(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%calving%calving_rate = inarray
  end subroutine glide_set_calving_rate

  subroutine glide_get_calving_thck(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%calving%calving_thck)
  end subroutine glide_get_calving_thck

  subroutine glide_set_calving_thck(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%calving%calving_thck = inarray/(thk0)
  end subroutine glide_set_calving_thck

  subroutine glide_get_cell_area(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (len0*len0)*(data%geometry%cell_area)
  end subroutine glide_get_cell_area

  subroutine glide_set_cell_area(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%cell_area = inarray/(len0*len0)
  end subroutine glide_set_cell_area

  subroutine glide_get_coulomb_c(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_physics%coulomb_c
  end subroutine glide_get_coulomb_c

  subroutine glide_set_coulomb_c(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%coulomb_c = inarray
  end subroutine glide_set_coulomb_c

  subroutine glide_get_coulomb_c_relax(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_physics%coulomb_c_relax
  end subroutine glide_get_coulomb_c_relax

  subroutine glide_set_coulomb_c_relax(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%coulomb_c_relax = inarray
  end subroutine glide_set_coulomb_c_relax

  subroutine glide_get_deltaT_ocn(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%ocean_data%deltaT_ocn
  end subroutine glide_get_deltaT_ocn


  subroutine glide_set_deltaT_ocn(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%ocean_data%deltaT_ocn = inarray
  end subroutine glide_set_deltaT_ocn


  subroutine glide_get_diff_cfl_dt(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%numerics%diff_cfl_dt
  end subroutine glide_get_diff_cfl_dt

  subroutine glide_set_diff_cfl_dt(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%numerics%diff_cfl_dt = inarray
  end subroutine glide_set_diff_cfl_dt

  subroutine glide_get_diffu(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_diffu)*(data%velocity%diffu)
  end subroutine glide_get_diffu

  subroutine glide_set_diffu(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%diffu = inarray/(scale_diffu)
  end subroutine glide_set_diffu

  subroutine glide_get_divu(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scyr)*(data%velocity%divu)
  end subroutine glide_get_divu

  subroutine glide_set_divu(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%divu = inarray/(scyr)
  end subroutine glide_set_divu

  subroutine glide_get_dthck_dt(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scyr)*(data%geometry%dthck_dt)
  end subroutine glide_get_dthck_dt

  subroutine glide_set_dthck_dt(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%dthck_dt = inarray/(scyr)
  end subroutine glide_set_dthck_dt

  subroutine glide_get_dthck_dt_obs(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%dthck_dt_obs
  end subroutine glide_get_dthck_dt_obs

  subroutine glide_set_dthck_dt_obs(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%dthck_dt_obs = inarray
  end subroutine glide_set_dthck_dt_obs

  subroutine glide_get_dthck_dt_obs_basin(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%dthck_dt_obs_basin
  end subroutine glide_get_dthck_dt_obs_basin

  subroutine glide_set_dthck_dt_obs_basin(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%dthck_dt_obs_basin = inarray
  end subroutine glide_set_dthck_dt_obs_basin

  subroutine glide_get_dthckdtm(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%geomderv%dthckdtm)
  end subroutine glide_get_dthckdtm

  subroutine glide_set_dthckdtm(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geomderv%dthckdtm = inarray/(scale_acab)
  end subroutine glide_set_dthckdtm

  subroutine glide_get_dusrfdtm(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_acab)*(data%geomderv%dusrfdtm)
  end subroutine glide_get_dusrfdtm

  subroutine glide_set_dusrfdtm(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geomderv%dusrfdtm = inarray/(scale_acab)
  end subroutine glide_set_dusrfdtm

  subroutine glide_get_effecpress(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_physics%effecpress
  end subroutine glide_get_effecpress

  subroutine glide_set_effecpress(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%effecpress = inarray
  end subroutine glide_set_effecpress

  subroutine glide_get_eps_eigen1(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scyr)*(data%calving%eps_eigen1)
  end subroutine glide_get_eps_eigen1

  subroutine glide_set_eps_eigen1(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%calving%eps_eigen1 = inarray/(scyr)
  end subroutine glide_set_eps_eigen1

  subroutine glide_get_eps_eigen2(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scyr)*(data%calving%eps_eigen2)
  end subroutine glide_get_eps_eigen2

  subroutine glide_set_eps_eigen2(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%calving%eps_eigen2 = inarray/(scyr)
  end subroutine glide_set_eps_eigen2

  subroutine glide_get_eus(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = (thk0)*(data%climate%eus)
  end subroutine glide_get_eus

  subroutine glide_set_eus(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%climate%eus = inarray/(thk0)
  end subroutine glide_set_eus

  subroutine glide_get_f_effecpress_bwat(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_physics%f_effecpress_bwat
  end subroutine glide_get_f_effecpress_bwat

  subroutine glide_set_f_effecpress_bwat(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%f_effecpress_bwat = inarray
  end subroutine glide_set_f_effecpress_bwat

  subroutine glide_get_f_effecpress_bwat_target(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_physics%f_effecpress_bwat_target
  end subroutine glide_get_f_effecpress_bwat_target

  subroutine glide_set_f_effecpress_bwat_target(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%f_effecpress_bwat_target = inarray
  end subroutine glide_set_f_effecpress_bwat_target

  subroutine glide_get_f_effecpress_bwatflx(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_physics%f_effecpress_bwatflx
  end subroutine glide_get_f_effecpress_bwatflx

  subroutine glide_set_f_effecpress_bwatflx(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%f_effecpress_bwatflx = inarray
  end subroutine glide_set_f_effecpress_bwatflx

  subroutine glide_get_f_effecpress_ocean_p(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_physics%f_effecpress_ocean_p
  end subroutine glide_get_f_effecpress_ocean_p

  subroutine glide_set_f_effecpress_ocean_p(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%f_effecpress_ocean_p = inarray
  end subroutine glide_set_f_effecpress_ocean_p

  subroutine glide_get_f_flotation(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%f_flotation
  end subroutine glide_get_f_flotation

  subroutine glide_set_f_flotation(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%f_flotation = inarray
  end subroutine glide_set_f_flotation

  subroutine glide_get_f_ground(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%f_ground
  end subroutine glide_get_f_ground

  subroutine glide_set_f_ground(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%f_ground = inarray
  end subroutine glide_set_f_ground

  subroutine glide_get_f_ground_cell(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%f_ground_cell
  end subroutine glide_get_f_ground_cell

  subroutine glide_set_f_ground_cell(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%f_ground_cell = inarray
  end subroutine glide_set_f_ground_cell

  subroutine glide_get_f_ground_obs(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%f_ground_obs
  end subroutine glide_get_f_ground_obs

  subroutine glide_set_f_ground_obs(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%f_ground_obs = inarray
  end subroutine glide_set_f_ground_obs

  subroutine glide_get_ff_invert_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%inversion%ff_invert_mask
  end subroutine glide_get_ff_invert_mask

  subroutine glide_set_ff_invert_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%inversion%ff_invert_mask = inarray
  end subroutine glide_set_ff_invert_mask

  subroutine glide_get_floating_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%floating_mask
  end subroutine glide_get_floating_mask

  subroutine glide_set_floating_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%geometry%floating_mask = inarray
  end subroutine glide_set_floating_mask

  subroutine glide_get_floating_thck_target(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%inversion%floating_thck_target)
  end subroutine glide_get_floating_thck_target

  subroutine glide_set_floating_thck_target(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%inversion%floating_thck_target = inarray/(thk0)
  end subroutine glide_set_floating_thck_target

  subroutine glide_get_flow_enhancement_factor(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%temper%flow_enhancement_factor
  end subroutine glide_get_flow_enhancement_factor

  subroutine glide_set_flow_enhancement_factor(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%temper%flow_enhancement_factor = inarray
  end subroutine glide_set_flow_enhancement_factor

  subroutine glide_get_gl_flux(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%gl_flux
  end subroutine glide_get_gl_flux

  subroutine glide_set_gl_flux(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%gl_flux = inarray
  end subroutine glide_set_gl_flux

  subroutine glide_get_gl_flux_east(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%gl_flux_east
  end subroutine glide_get_gl_flux_east

  subroutine glide_set_gl_flux_east(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%gl_flux_east = inarray
  end subroutine glide_set_gl_flux_east

  subroutine glide_get_gl_flux_north(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%gl_flux_north
  end subroutine glide_get_gl_flux_north

  subroutine glide_set_gl_flux_north(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%gl_flux_north = inarray
  end subroutine glide_set_gl_flux_north

  subroutine glide_get_gravity(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = grav
  end subroutine glide_get_gravity

  subroutine glide_set_gravity(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_gravity

  subroutine glide_get_grounded_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%grounded_mask
  end subroutine glide_get_grounded_mask

  subroutine glide_set_grounded_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%geometry%grounded_mask = inarray
  end subroutine glide_set_grounded_mask

  subroutine glide_get_head(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1)*(data%basal_hydro%head)
  end subroutine glide_get_head

  subroutine glide_set_head(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_hydro%head = inarray/(1)
  end subroutine glide_set_head

  subroutine glide_get_iarea(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%geometry%iarea
  end subroutine glide_get_iarea

  subroutine glide_set_iarea(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%iarea = inarray
  end subroutine glide_set_iarea

  subroutine glide_get_iareaf(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%geometry%iareaf
  end subroutine glide_get_iareaf

  subroutine glide_set_iareaf(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%iareaf = inarray
  end subroutine glide_set_iareaf

  subroutine glide_get_iareag(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%geometry%iareag
  end subroutine glide_get_iareag

  subroutine glide_set_iareag(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%iareag = inarray
  end subroutine glide_set_iareag

  subroutine glide_get_ice_cap_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%ice_cap_mask
  end subroutine glide_get_ice_cap_mask

  subroutine glide_set_ice_cap_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%geometry%ice_cap_mask = inarray
  end subroutine glide_set_ice_cap_mask

  subroutine glide_get_ice_domain_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%general%ice_domain_mask
  end subroutine glide_get_ice_domain_mask

  subroutine glide_set_ice_domain_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%general%ice_domain_mask = inarray
  end subroutine glide_set_ice_domain_mask

  subroutine glide_get_ice_fraction_retreat_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%ice_fraction_retreat_mask
  end subroutine glide_get_ice_fraction_retreat_mask

  subroutine glide_set_ice_fraction_retreat_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%ice_fraction_retreat_mask = inarray
  end subroutine glide_set_ice_fraction_retreat_mask

  subroutine glide_get_ice_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%ice_mask
  end subroutine glide_get_ice_mask

  subroutine glide_set_ice_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%geometry%ice_mask = inarray
  end subroutine glide_set_ice_mask

  subroutine glide_get_ice_mask_stag(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%ice_mask_stag
  end subroutine glide_get_ice_mask_stag

  subroutine glide_set_ice_mask_stag(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%geometry%ice_mask_stag = inarray
  end subroutine glide_set_ice_mask_stag

  subroutine glide_get_ice_sheet_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%ice_sheet_mask
  end subroutine glide_get_ice_sheet_mask

  subroutine glide_set_ice_sheet_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%geometry%ice_sheet_mask = inarray
  end subroutine glide_set_ice_sheet_mask

  subroutine glide_get_ice_specific_heat(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = shci
  end subroutine glide_get_ice_specific_heat

  subroutine glide_set_ice_specific_heat(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_ice_specific_heat

  subroutine glide_get_ice_thermal_conductivity(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = coni
  end subroutine glide_get_ice_thermal_conductivity

  subroutine glide_set_ice_thermal_conductivity(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_ice_thermal_conductivity

  subroutine glide_get_imass(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%geometry%imass
  end subroutine glide_get_imass

  subroutine glide_set_imass(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%imass = inarray
  end subroutine glide_set_imass

  subroutine glide_get_imass_above_flotation(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%geometry%imass_above_flotation
  end subroutine glide_get_imass_above_flotation

  subroutine glide_set_imass_above_flotation(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%imass_above_flotation = inarray
  end subroutine glide_set_imass_above_flotation

  subroutine glide_get_ivol(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%geometry%ivol
  end subroutine glide_get_ivol

  subroutine glide_set_ivol(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%ivol = inarray
  end subroutine glide_set_ivol

  subroutine glide_get_kinbcmask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%velocity%kinbcmask(:,:)
  end subroutine glide_get_kinbcmask

  subroutine glide_set_kinbcmask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%velocity%kinbcmask(:,:) = inarray
  end subroutine glide_set_kinbcmask

  subroutine glide_get_lat(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%general%lat
  end subroutine glide_get_lat

  subroutine glide_set_lat(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%general%lat = inarray
  end subroutine glide_set_lat

  subroutine glide_get_load(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%isostasy%load)
  end subroutine glide_get_load

  subroutine glide_set_load(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%isostasy%load = inarray/(thk0)
  end subroutine glide_set_load

  subroutine glide_get_lon(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%general%lon
  end subroutine glide_get_lon

  subroutine glide_set_lon(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%general%lon = inarray
  end subroutine glide_set_lon

  subroutine glide_get_lsurf(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geometry%lsrf)
  end subroutine glide_get_lsurf

  subroutine glide_set_lsurf(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%lsrf = inarray/(thk0)
  end subroutine glide_set_lsurf

  subroutine glide_get_marine_connection_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%marine_connection_mask
  end subroutine glide_get_marine_connection_mask

  subroutine glide_set_marine_connection_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%geometry%marine_connection_mask = inarray
  end subroutine glide_set_marine_connection_mask

  subroutine glide_get_marine_connection_mask_isolated(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%marine_connection_mask_isolated
  end subroutine glide_get_marine_connection_mask_isolated

  subroutine glide_set_marine_connection_mask_isolated(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%geometry%marine_connection_mask_isolated = inarray
  end subroutine glide_set_marine_connection_mask_isolated

  subroutine glide_get_overwrite_acab_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%climate%overwrite_acab_mask
  end subroutine glide_get_overwrite_acab_mask

  subroutine glide_set_overwrite_acab_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%climate%overwrite_acab_mask = inarray
  end subroutine glide_set_overwrite_acab_mask

  subroutine glide_get_powerlaw_c(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%basal_physics%powerlaw_c
  end subroutine glide_get_powerlaw_c

  subroutine glide_set_powerlaw_c(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%powerlaw_c = inarray
  end subroutine glide_set_powerlaw_c

  subroutine glide_get_reference_thck(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geometry%reference_thck)
  end subroutine glide_get_reference_thck

  subroutine glide_set_reference_thck(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%reference_thck = inarray/(thk0)
  end subroutine glide_set_reference_thck

  subroutine glide_get_relx(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%isostasy%relx)
  end subroutine glide_get_relx

  subroutine glide_set_relx(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%isostasy%relx = inarray/(thk0)
  end subroutine glide_set_relx

  subroutine glide_get_rho_ice(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = rhoi
  end subroutine glide_get_rho_ice

  subroutine glide_set_rho_ice(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_rho_ice

  subroutine glide_get_rho_seawater(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = rhoo
  end subroutine glide_get_rho_seawater

  subroutine glide_set_rho_seawater(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_rho_seawater

  subroutine glide_get_seconds_per_year(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = scyr
  end subroutine glide_get_seconds_per_year

  subroutine glide_set_seconds_per_year(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

!  no rescaling here
  end subroutine glide_set_seconds_per_year

  subroutine glide_get_sfc_mbal_flux(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%sfc_mbal_flux
  end subroutine glide_get_sfc_mbal_flux

  subroutine glide_set_sfc_mbal_flux(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%sfc_mbal_flux = inarray
  end subroutine glide_set_sfc_mbal_flux

  subroutine glide_get_smb(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1.0)*(data%climate%smb)
  end subroutine glide_get_smb

  subroutine glide_set_smb(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%smb = inarray/(1.0)
  end subroutine glide_set_smb

  subroutine glide_get_smb_anomaly(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1.0)*(data%climate%smb_anomaly)
  end subroutine glide_get_smb_anomaly

  subroutine glide_set_smb_anomaly(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%smb_anomaly = inarray/(1.0)
  end subroutine glide_set_smb_anomaly

  subroutine glide_get_smb_gradz(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1.0/thk0)*(data%climate%smb_gradz)
  end subroutine glide_get_smb_gradz

  subroutine glide_set_smb_gradz(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%smb_gradz = inarray/(1.0/thk0)
  end subroutine glide_set_smb_gradz

  subroutine glide_get_smb_levels(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:), intent(out) :: outarray

    outarray = (thk0)*(data%climate%smb_levels)
  end subroutine glide_get_smb_levels

  subroutine glide_set_smb_levels(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:), intent(in) :: inarray

    data%climate%smb_levels = inarray/(thk0)
  end subroutine glide_set_smb_levels

  subroutine glide_get_smb_ref(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1.0)*(data%climate%smb_ref)
  end subroutine glide_get_smb_ref

  subroutine glide_set_smb_ref(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%smb_ref = inarray/(1.0)
  end subroutine glide_set_smb_ref

  subroutine glide_get_smb_reference_usrf(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%climate%smb_reference_usrf)
  end subroutine glide_get_smb_reference_usrf

  subroutine glide_set_smb_reference_usrf(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%climate%smb_reference_usrf = inarray/(thk0)
  end subroutine glide_set_smb_reference_usrf

  subroutine glide_get_soft(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_btrc)*(data%velocity%bed_softness)
  end subroutine glide_get_soft

  subroutine glide_set_soft(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%bed_softness = inarray/(scale_btrc)
  end subroutine glide_set_soft

  subroutine glide_get_stagthk(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geomderv%stagthck)
  end subroutine glide_get_stagthk

  subroutine glide_set_stagthk(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geomderv%stagthck = inarray/(thk0)
  end subroutine glide_set_stagthk

  subroutine glide_get_tau_c(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1e-3)*(data%basal_physics%tau_c)
  end subroutine glide_get_tau_c

  subroutine glide_set_tau_c(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%tau_c = inarray/(1e-3)
  end subroutine glide_set_tau_c

  subroutine glide_get_tau_eff_calving(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%calving%tau_eff
  end subroutine glide_get_tau_eff_calving

  subroutine glide_set_tau_eff_calving(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%calving%tau_eff = inarray
  end subroutine glide_set_tau_eff_calving

  subroutine glide_get_tau_eigen1(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%calving%tau_eigen1
  end subroutine glide_get_tau_eigen1

  subroutine glide_set_tau_eigen1(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%calving%tau_eigen1 = inarray
  end subroutine glide_set_tau_eigen1

  subroutine glide_get_tau_eigen2(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%calving%tau_eigen2
  end subroutine glide_get_tau_eigen2

  subroutine glide_set_tau_eigen2(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%calving%tau_eigen2 = inarray
  end subroutine glide_set_tau_eigen2

  subroutine glide_get_taudx(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%taudx(:,:))
  end subroutine glide_get_taudx

  subroutine glide_set_taudx(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%taudx(:,:) = inarray/(scale_tau)
  end subroutine glide_set_taudx

  subroutine glide_get_taudy(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%stress%taudy(:,:))
  end subroutine glide_get_taudy

  subroutine glide_set_taudy(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%stress%taudy(:,:) = inarray/(scale_tau)
  end subroutine glide_set_taudy

  subroutine glide_get_tauf(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_tau)*(data%basal_physics%mintauf)
  end subroutine glide_get_tauf

  subroutine glide_set_tauf(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%basal_physics%mintauf = inarray/(scale_tau)
  end subroutine glide_set_tauf

  subroutine glide_get_taux(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1e-3*thk0*thk0/len0)*(data%velocity%tau_x)
  end subroutine glide_get_taux

  subroutine glide_set_taux(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%tau_x = inarray/(1e-3*thk0*thk0/len0)
  end subroutine glide_set_taux

  subroutine glide_get_tauy(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (1e-3*thk0*thk0/len0)*(data%velocity%tau_y)
  end subroutine glide_get_tauy

  subroutine glide_set_tauy(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%tau_y = inarray/(1e-3*thk0*thk0/len0)
  end subroutine glide_set_tauy

  subroutine glide_get_thck_calving_threshold(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%calving%thck_calving_threshold
  end subroutine glide_get_thck_calving_threshold

  subroutine glide_set_thck_calving_threshold(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%calving%thck_calving_threshold = inarray
  end subroutine glide_set_thck_calving_threshold

  subroutine glide_get_thermal_forcing_lsrf(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = data%ocean_data%thermal_forcing_lsrf(:,:)
  end subroutine glide_get_thermal_forcing_lsrf

  subroutine glide_get_thk(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geometry%thck)
  end subroutine glide_get_thk

  subroutine glide_set_thk(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%thck = inarray/(thk0)
  end subroutine glide_set_thk

  subroutine glide_get_thkmask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%geometry%thkmask
  end subroutine glide_get_thkmask

  subroutine glide_set_thkmask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%geometry%thkmask = inarray
  end subroutine glide_set_thkmask

  subroutine glide_get_topg(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geometry%topg)
  end subroutine glide_get_topg

  subroutine glide_set_topg(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%topg = inarray/(thk0)
  end subroutine glide_set_topg

  subroutine glide_get_topg_stdev(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geometry%topg_stdev)
  end subroutine glide_get_topg_stdev

  subroutine glide_set_topg_stdev(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%topg_stdev = inarray/(thk0)
  end subroutine glide_set_topg_stdev

  subroutine glide_get_total_bmb_flux(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%geometry%total_bmb_flux
  end subroutine glide_get_total_bmb_flux

  subroutine glide_set_total_bmb_flux(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%total_bmb_flux = inarray
  end subroutine glide_set_total_bmb_flux

  subroutine glide_get_total_calving_flux(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%geometry%total_calving_flux
  end subroutine glide_get_total_calving_flux

  subroutine glide_set_total_calving_flux(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%total_calving_flux = inarray
  end subroutine glide_set_total_calving_flux

  subroutine glide_get_total_gl_flux(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%geometry%total_gl_flux
  end subroutine glide_get_total_gl_flux

  subroutine glide_set_total_gl_flux(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%total_gl_flux = inarray
  end subroutine glide_set_total_gl_flux

  subroutine glide_get_total_smb_flux(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(out) :: outarray

    outarray = data%geometry%total_smb_flux
  end subroutine glide_get_total_smb_flux

  subroutine glide_set_total_smb_flux(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), intent(in) :: inarray

    data%geometry%total_smb_flux = inarray
  end subroutine glide_set_total_smb_flux

  subroutine glide_get_ubas(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%uvel(data%general%upn,1:data%general%ewn-1,1:data%general%nsn-1))
  end subroutine glide_get_ubas

  subroutine glide_get_uflx(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uflx)*(data%velocity%uflx)
  end subroutine glide_get_uflx

  subroutine glide_set_uflx(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%uflx = inarray/(scale_uflx)
  end subroutine glide_set_uflx

  subroutine glide_get_unstagbeta(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_beta)*(data%velocity%unstagbeta)
  end subroutine glide_get_unstagbeta

  subroutine glide_set_unstagbeta(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%unstagbeta = inarray/(scale_beta)
  end subroutine glide_set_unstagbeta

  subroutine glide_get_usfc(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%uvel(1,1:data%general%ewn-1,1:data%general%nsn-1))
  end subroutine glide_get_usfc

  subroutine glide_get_usfc_obs(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%usfc_obs)
  end subroutine glide_get_usfc_obs

  subroutine glide_set_usfc_obs(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%usfc_obs = inarray/(scale_uvel)
  end subroutine glide_set_usfc_obs

  subroutine glide_get_usrf_obs(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geometry%usrf_obs)
  end subroutine glide_get_usrf_obs

  subroutine glide_set_usrf_obs(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%usrf_obs = inarray/(thk0)
  end subroutine glide_set_usrf_obs

  subroutine glide_get_usurf(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (thk0)*(data%geometry%usrf)
  end subroutine glide_get_usurf

  subroutine glide_set_usurf(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%geometry%usrf = inarray/(thk0)
  end subroutine glide_set_usurf

  subroutine glide_get_uvel_2d(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%uvel_2d(:,:))
  end subroutine glide_get_uvel_2d

  subroutine glide_set_uvel_2d(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%uvel_2d(:,:) = inarray/(scale_uvel)
  end subroutine glide_set_uvel_2d

  subroutine glide_get_uvel_2d_extend(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%uvel_2d_extend(:,:))
  end subroutine glide_get_uvel_2d_extend

  subroutine glide_set_uvel_2d_extend(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%uvel_2d_extend(:,:) = inarray/(scale_uvel)
  end subroutine glide_set_uvel_2d_extend

  subroutine glide_get_uvel_mean(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%uvel_mean)
  end subroutine glide_get_uvel_mean

  subroutine glide_set_uvel_mean(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%uvel_mean = inarray/(scale_uvel)
  end subroutine glide_set_uvel_mean

  subroutine glide_get_vbas(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%vvel(data%general%upn,1:data%general%ewn-1,1:data%general%nsn-1))
  end subroutine glide_get_vbas

  subroutine glide_get_velo_sfc_obs(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%velo_sfc_obs)
  end subroutine glide_get_velo_sfc_obs

  subroutine glide_set_velo_sfc_obs(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%velo_sfc_obs = inarray/(scale_uvel)
  end subroutine glide_set_velo_sfc_obs

  subroutine glide_get_vflx(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uflx)*(data%velocity%vflx)
  end subroutine glide_get_vflx

  subroutine glide_set_vflx(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%vflx = inarray/(scale_uflx)
  end subroutine glide_set_vflx

  subroutine glide_get_vsfc(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%vvel(1,1:data%general%ewn-1,1:data%general%nsn-1))
  end subroutine glide_get_vsfc

  subroutine glide_get_vsfc_obs(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%vsfc_obs)
  end subroutine glide_get_vsfc_obs

  subroutine glide_set_vsfc_obs(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%vsfc_obs = inarray/(scale_uvel)
  end subroutine glide_set_vsfc_obs

  subroutine glide_get_vvel_2d(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%vvel_2d(:,:))
  end subroutine glide_get_vvel_2d

  subroutine glide_set_vvel_2d(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%vvel_2d(:,:) = inarray/(scale_uvel)
  end subroutine glide_set_vvel_2d

  subroutine glide_get_vvel_2d_extend(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%vvel_2d_extend(:,:))
  end subroutine glide_get_vvel_2d_extend

  subroutine glide_set_vvel_2d_extend(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%vvel_2d_extend(:,:) = inarray/(scale_uvel)
  end subroutine glide_set_vvel_2d_extend

  subroutine glide_get_vvel_mean(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%vvel_mean)
  end subroutine glide_get_vvel_mean

  subroutine glide_set_vvel_mean(data,inarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(in) :: inarray

    data%velocity%vvel_mean = inarray/(scale_uvel)
  end subroutine glide_set_vvel_mean

  subroutine glide_get_warm_ocean_mask(data,outarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(out) :: outarray

    outarray = data%basal_melt%warm_ocean_mask
  end subroutine glide_get_warm_ocean_mask

  subroutine glide_set_warm_ocean_mask(data,inarray)
    implicit none
    type(glide_global_type) :: data
    integer, dimension(:,:), intent(in) :: inarray

    data%basal_melt%warm_ocean_mask = inarray
  end subroutine glide_set_warm_ocean_mask

  subroutine glide_get_wbas(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%wvel(data%general%upn,1:data%general%ewn,1:data%general%nsn))
  end subroutine glide_get_wbas

  subroutine glide_get_wsfc(data,outarray)
    implicit none
    type(glide_global_type) :: data
    real(dp), dimension(:,:), intent(out) :: outarray

    outarray = (scale_uvel)*(data%velocity%wvel(1,1:data%general%ewn,1:data%general%nsn))
  end subroutine glide_get_wsfc


end module glide_io
