!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_ncio.F90 - part of the Community Ice Sheet Model (CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#define NCO outfile%nc
#define NCI infile%nc

module glimmer_ncio
  !> module for common netCDF I/O
  !> written by Magnus Hagdorn, 2004

  use glimmer_ncdf
  use cism_parallel, only: parallel_type, parallel_create, parallel_open, parallel_put_var, parallel_get_var, &
       parallel_put_att, parallel_def_var, parallel_def_dim, parallel_inq_varid, parallel_inq_dimid,  &
       parallel_inquire_dimension, parallel_redef, parallel_enddef, parallel_sync, parallel_get_att

  implicit none

  ! All routines in this module are public

  integer,parameter,private :: msglen=512
  integer,parameter,private :: tavg_list_len=4096   ! length of the list of time-average stream names
  
  interface glimmer_nc_get_var
     module procedure glimmer_nc_get_var_integer_2d
     module procedure glimmer_nc_get_var_real8_2d
  end interface

  logical, parameter :: verbose_ncio = .false.

contains

  !*****************************************************************************
  ! netCDF output
  !*****************************************************************************  
  subroutine openall_out(model,outfiles)

    !> open all netCDF files for output
    use glide_types
    use glimmer_ncdf
    use glimmer_filenames, only: process_path

    implicit none

    type(glide_global_type) :: model
    type(glimmer_nc_output), pointer, optional :: outfiles
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc
    integer :: status

    if (present(outfiles)) then
       oc => outfiles
    else
       oc => model%funits%out_first
    end if

    do while(associated(oc))

       if (oc%one_file_per_write) then

          ! No file at initialization; a new file is created at each write

       elseif (oc%append) then   ! assume the file exists, and reopen it

          call glimmer_nc_openappend(oc,model)

       elseif (model%options%is_restart == STANDARD_RESTART) then   ! reopen the file if it exists

          status = parallel_open(process_path(oc%nc%filename),NF90_WRITE,oc%nc%id)

          if (status == NF90_NOERR) then  ! file exists and is now open; append it

             oc%append = .true.
             call glimmer_nc_openappend(oc, model, already_open_in=.true.)

          else  ! file does not exist; create it

             call glimmer_nc_createfile(oc, model)

          endif

       else  ! assume the file does not exist; create it
             ! Note: For hybrid restarts, the file is created at initialization

          call glimmer_nc_createfile(oc, model)

       end if

       ! Set the start of the first averaging interval (used for tavg time bounds)
       call glimmer_nc_checkwrite_init(oc, model%numerics%tstart, &
                                       tstep_count = model%numerics%tstep_count)

       oc => oc%next

    end do

  end subroutine openall_out

  !------------------------------------------------------------------------------

  subroutine closeall_out(model,outfiles)

    !> close all netCDF files for output
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    type(glimmer_nc_output),pointer,optional :: outfiles

    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc => model%funits%out_first
    end if

    do while(associated(oc))
       oc => delete(oc)
    end do
    if (.not.present(outfiles)) model%funits%out_first=>NULL()

  end subroutine closeall_out

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_openappend(outfile, model, &
                                   already_open_in)

    !> open netCDF file for appending
    use glimmer_log
    use glide_types
    use glimmer_map_CFproj
    use glimmer_map_types
    use glimmer_filenames

    implicit none

    type(glimmer_nc_output), pointer :: outfile       !> structure containing output netCDF descriptor
    type(glide_global_type) :: model                  !> the model instance
    logical, intent(in), optional :: already_open_in  !> if true, then the file is already open 

    ! local variables
    integer :: status, timedimid, ntime
    integer :: nstreams                     ! number of time-average output streams
    character(len=tavg_list_len) :: stream_list, file_stream_list
    character(len=msglen) :: message
    logical :: already_open   ! if true, the file is already open

    if (present(already_open_in)) then
       already_open = already_open_in
    else
       already_open = .false.
    endif

    ! open the existing netCDF file, if not already open
    if (.not. already_open) then
       status = parallel_open(process_path(NCO%filename),NF90_WRITE,NCO%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
    endif
    NCO%file_open = .true.

    call write_log_div
    write(message,*) 'Reopening file ',trim(process_path(NCO%filename)),' for output; '
    call write_log(trim(message))

    ! Find out when last time slice was written
    status = parallel_inq_dimid(NCO%id,'time',timedimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inquire_dimension(NCO%id,timedimid,len=ntime)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! Set timecounter
    outfile%timecounter = ntime+1

    write(message,*) '  Write every ', outfile%freq, ' years'
    call write_log(trim(message))

    ! Get time-related varids
    status = parallel_inq_varid(NCO%id,glimmer_nc_internal_time_varname,NCO%internal_timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inq_varid(NCO%id,glimmer_nc_time_varname,NCO%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inq_varid(NCO%id,glimmer_nc_tstep_count_varname,NCO%tstep_count_var)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! For time-average files, get the time bounds varids
    ! Note: This test for '_tavg' matches the test in glimmer_nc_createfile.
    NCO%time_bounds = .false.
    if (index(NCO%vars,'_tavg') /= 0) then
       status = parallel_inq_varid(NCO%id,glimmer_nc_internal_timebounds_varname,NCO%internal_timebounds_var)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_inq_varid(NCO%id,glimmer_nc_timebounds_varname,NCO%timebounds_var)
       call nc_errorhandle(__FILE__,__LINE__,status)
       NCO%time_bounds = .true.
    end if

    ! For a restart file, find the variables that hold the averaging state of the time-average
    ! output streams (see glimmer_nc_createfile). The averaging state is written to this file only
    ! if the file lists the same streams as the current run.
    NCO%tavg_state = .false.
    status = parallel_inq_dimid(NCO%id, 'tavgstream', NCO%tavgstream_dim)
    if (status == NF90_NOERR) then
       call glimmer_nc_tavg_streams(model, nstreams, stream_list)
       file_stream_list = ''
       status = parallel_get_att(NCO%id, NF90_GLOBAL, 'tavg_streams', file_stream_list)
       if (status == NF90_NOERR .and. trim(file_stream_list) == trim(stream_list)) then
          status = parallel_inq_varid(NCO%id, 'tavg_total_time', NCO%tavg_total_time_var)
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_inq_varid(NCO%id, 'tavg_start_time', NCO%tavg_start_time_var)
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_inq_varid(NCO%id, 'tavg_start_external_time', NCO%tavg_start_external_time_var)
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_inq_varid(NCO%id, 'tavg_accum_tstep_count', NCO%tavg_accum_tstep_count_var)
          call nc_errorhandle(__FILE__,__LINE__,status)
          NCO%tavg_state = .true.
       else
          call write_log('The time-average output streams listed in '//trim(process_path(NCO%filename))// &
               ' ('//trim(file_stream_list)//') differ from those of this run ('//trim(stream_list)// &
               '); the averaging state will not be written to this file', GM_WARNING)
       end if
    end if

    ! Put dataset into define mode
    status = parallel_redef(NCO%id)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! level and staglevel dimension
    NCO%nlevel = model%general%upn
    NCO%nstaglevel = model%general%upn-1
    NCO%nstagwbndlevel = model%general%upn ! MJH this is the max index, not the size

    ! Note: The following dimension lengths are set to 1 by default,
    !       but can be increased depending on the config options.

    ! vertical coordinate for ocean data
    NCO%nzocn = model%ocean_data%nzocn

    ! vertical coordinate for atmosphere data
    NCO%nzatm = model%climate%nzatm

    ! glacier ID coordinate for glacier data
    NCO%nglacier = model%glacier%nglacier

    ! basin coordinate for basin data
    NCO%nbasin = model%ocean_data%nbasin

    ! coordinate for calvingMIP axis data
    NCO%naxis = model%calving%naxis

  end subroutine glimmer_nc_openappend

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_createfile(outfile, model, external_baseline_year, external_time_units)

    !> create a new netCDF file
    use glimmer_log
    use glide_types
    use glimmer_map_CFproj
    use glimmer_map_types
    use glimmer_filenames
    implicit none

    type(glimmer_nc_output), pointer :: outfile     !> structure containing output netCDF descriptor
    type(glide_global_type) :: model                !> the model instance
    integer, intent(in), optional :: &
         external_baseline_year                     !> baseline year for external time; default = internal_baseline_year
    character(len=*), intent(in), optional :: &     !> external time units; default = internal_time_units
         external_time_units

    ! local variables
    integer, parameter :: time_units_len = 128
    integer status
    integer mapid
    integer :: pos
    integer :: sub_external_baseline_year  ! local version of external_baseline_year
    character(len=:), allocatable :: sub_external_time_units  ! local version of external_time_units
    character(len=4) :: year_str
    character(len=time_units_len) :: internal_time_units_str
    character(len=time_units_len) :: time_units_str
    integer :: nstreams                     ! number of time-average output streams
    character(len=tavg_list_len) :: stream_list
    character(len=msglen) message

    ! Note: The internal baseline year and units are hardcoded.
    !       The default baseline year used to be year 1, implying that t = 0 corresponds to Jan. 1 of year 1,
    !        hence t = 1950.0 is Jan. 1 of year 1951.
    !       With a baseline year of 0, t = 1950.0 corresponds to Jan. 1, 1950, which is more intuitive.

    integer, parameter :: internal_baseline_year = 0
    ! Note: 'common_years' (with an underscore) is the UDUNITS name for years of exactly 365 days.
    !       With a space ('common years'), cftime and xarray cannot decode the time variables.
    character(len=*), parameter :: internal_time_units = 'common_years'  ! common year = year of exactly 365 days

    if (present(external_baseline_year)) then
       sub_external_baseline_year = external_baseline_year
    else
       sub_external_baseline_year = internal_baseline_year
    end if

    if (present(external_time_units)) then
       sub_external_time_units = external_time_units
    else
       sub_external_time_units = internal_time_units
    endif

    ! create new netCDF file
    !WHL - Changed the following line to support large netCDF output files
!!    status = parallel_create(process_path(NCO%filename),NF90_CLOBBER,NCO%id)
    status = parallel_create(process_path(NCO%filename), ior(NF90_CLOBBER,NF90_64BIT_OFFSET), NCO%id)
    call nc_errorhandle(__FILE__,__LINE__,status)
    NCO%file_open = .true.
    call write_log_div
    write(message,*) 'Opening file ', trim(process_path(NCO%filename)), ' for output; '
    call write_log(trim(message))

    if (outfile%write_init) then
       write(message,*) '  Write output at start of run and every ', outfile%freq, ' years'
    else
       write(message,*) '  Write output every ', outfile%freq, ' years'
    endif
    call write_log(trim(message))

    if (outfile%end_write < glimmer_nc_max_time) then
       write(message,*) '  Stop writing at ', outfile%end_write
       call write_log(trim(message))
    end if
    NCO%define_mode=.TRUE.

    ! write meta data
    status = parallel_put_att(NCO%id, NF90_GLOBAL, 'Conventions', "CF-1.3")
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'title',trim(outfile%metadata%title))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'institution',trim(outfile%metadata%institution))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'source',trim(outfile%metadata%source))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'history',trim(outfile%metadata%history))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'references',trim(outfile%metadata%references))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'comment',trim(outfile%metadata%comment))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'configuration',trim(outfile%metadata%config))
    call nc_errorhandle(__FILE__,__LINE__,status)
  
    ! define time dimension
    status = parallel_def_dim(NCO%id,'time',NF90_UNLIMITED,NCO%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! See comments in glimmer_ncdf explaining the multiple time-related variables.
    call write_log('Creating variables internal_time, time, and tstep_count')

    ! define the internal_time variable

    status = parallel_def_var(NCO%id,glimmer_nc_internal_time_varname,&
         outfile%time_xtype,(/NCO%timedim/),NCO%internal_timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NCO%internal_timevar, 'long_name', 'internal time')

    ! Internal time units are hardcoded as common years (i.e., years of exactly 365 days).
    ! The baseline year is in YYYY format, with additional digits as needed for years > 9999.
    write(year_str,'(i0.4)') internal_baseline_year
    internal_time_units_str = internal_time_units // ' since ' // year_str // '-01-01 00:00:00'
    status = parallel_put_att(NCO%id, NCO%internal_timevar, 'units', internal_time_units_str)

    ! CISM currently assumes a noleap calendar - exactly 365 days. For now, we hard-code
    ! this assumption in the units (by hard-coding that we're using units of common_year:
    ! CF/Udunits defines common_year to be 365 days, whereas year means 365.242198781
    ! days) and the calendar attribute.
    status = parallel_put_att(NCO%id, NCO%internal_timevar, 'calendar', 'noleap')

    ! define the time variable
    ! By default, 'time' has the same properties as internal_time,
    ! but these can be overwritten by passing in an external time (e.g., from CESM)
    ! with different units and/or baseline_year.

    status = parallel_def_var(NCO%id,glimmer_nc_time_varname,&
         outfile%time_xtype,(/NCO%timedim/),NCO%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NCO%timevar, 'long_name', 'time')
    write(year_str,'(i0.4)') sub_external_baseline_year
    time_units_str = sub_external_time_units // ' since ' // year_str // '-01-01 00:00:00'
    status = parallel_put_att(NCO%id, NCO%timevar, 'units', time_units_str)
    status = parallel_put_att(NCO%id, NCO%timevar, 'calendar', 'noleap')

    ! define the tstep_count variable

    status = parallel_def_var(NCO%id,glimmer_nc_tstep_count_varname,&
         NF90_INT,(/NCO%timedim/),NCO%tstep_count_var)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NCO%tstep_count_var, 'long_name', &
         'Time step count')
    status = parallel_put_att(NCO%id, NCO%tstep_count_var, 'units', '-')

    ! If this is a time-average file, then add metadata for time bounds
    ! time_bounds has dimension (time,2) since there are two values (start and end) per time slice.

    ! Note: This test is done before NAME_io_create expands the 'restart' keyword.
    !       Restart files do not get time bounds, even if the expanded restart variable list
    !       includes variables with the '_tavg' suffix (e.g., the GLAD coupling fluxes rofi_tavg).

    NCO%time_bounds = .false.
    pos = index(NCO%vars,"_tavg")
    if (pos.ne.0) then  ! this is a time-average file

       NCO%time_bounds = .true.

       if (verbose_ncio .and. main_task) &
            write(iulog,*) 'Create time_bounds for file ', trim(NCO%filename)

       ! define a time bounds dimension ('nbnd', following CESM convention)
       status = parallel_def_dim(NCO%id,'nbnd',2,NCO%nbnd_dim)
       call nc_errorhandle(__FILE__,__LINE__,status)

       ! internal time bounds
       status = parallel_def_var(NCO%id,glimmer_nc_internal_timebounds_varname,&
            outfile%time_xtype,(/NCO%nbnd_dim,NCO%timedim/),NCO%internal_timebounds_var)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, NCO%internal_timevar, 'bounds', glimmer_nc_internal_timebounds_varname)
       status = parallel_put_att(NCO%id, NCO%internal_timebounds_var, 'long_name', &
            'internal time interval endpoints')
       status = parallel_put_att(NCO%id, NCO%internal_timebounds_var, 'units', &
            internal_time_units_str)
       status = parallel_put_att(NCO%id, NCO%internal_timebounds_var, 'calendar', 'noleap')

       ! external time bounds
       status = parallel_def_var(NCO%id,glimmer_nc_timebounds_varname,&
         outfile%time_xtype,(/NCO%nbnd_dim,NCO%timedim/),NCO%timebounds_var)
       call nc_errorhandle(__FILE__,__LINE__,status)
       status = parallel_put_att(NCO%id, NCO%timevar, 'bounds', glimmer_nc_timebounds_varname)
       status = parallel_put_att(NCO%id, NCO%timebounds_var, 'long_name', &
            'time interval endpoints')
       status = parallel_put_att(NCO%id, NCO%timebounds_var, 'units', time_units_str)
       status = parallel_put_att(NCO%id, NCO%timebounds_var, 'calendar', 'noleap')

    end if

    ! If this is a restart file, define variables to hold the averaging state of each time-average
    ! output stream, so that the averages continue exactly across a standard restart (restart = 1).
    ! (The running sums themselves are written as '_tavg_sum' variables; see NAME_io_createall.)
    ! The streams are identified by name, in order, in the global attribute 'tavg_streams'.
    ! Note: As for the time bounds above, this test is done before the 'restart' keyword is expanded.
    NCO%tavg_state = .false.
    if (glimmer_nc_is_restart(outfile)) then
       call glimmer_nc_tavg_streams(model, nstreams, stream_list)
       if (nstreams > 0) then
          NCO%tavg_state = .true.
          status = parallel_def_dim(NCO%id, 'tavgstream', nstreams, NCO%tavgstream_dim)
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_put_att(NCO%id, NF90_GLOBAL, 'tavg_streams', trim(stream_list))
          call nc_errorhandle(__FILE__,__LINE__,status)

          status = parallel_def_var(NCO%id, 'tavg_total_time', NF90_DOUBLE, &
               (/NCO%tavgstream_dim, NCO%timedim/), NCO%tavg_total_time_var)
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_put_att(NCO%id, NCO%tavg_total_time_var, 'long_name', &
               'time accumulated in the current averaging interval, for each time-average stream')
          status = parallel_put_att(NCO%id, NCO%tavg_total_time_var, 'units', trim(internal_time_units))

          status = parallel_def_var(NCO%id, 'tavg_start_time', NF90_DOUBLE, &
               (/NCO%tavgstream_dim, NCO%timedim/), NCO%tavg_start_time_var)
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_put_att(NCO%id, NCO%tavg_start_time_var, 'long_name', &
               'start of the current averaging interval (internal time), for each time-average stream')
          status = parallel_put_att(NCO%id, NCO%tavg_start_time_var, 'units', trim(internal_time_units_str))

          status = parallel_def_var(NCO%id, 'tavg_start_external_time', NF90_DOUBLE, &
               (/NCO%tavgstream_dim, NCO%timedim/), NCO%tavg_start_external_time_var)
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_put_att(NCO%id, NCO%tavg_start_external_time_var, 'long_name', &
               'start of the current averaging interval (external time), for each time-average stream')

          status = parallel_def_var(NCO%id, 'tavg_accum_tstep_count', NF90_INT, &
               (/NCO%tavgstream_dim, NCO%timedim/), NCO%tavg_accum_tstep_count_var)
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_put_att(NCO%id, NCO%tavg_accum_tstep_count_var, 'long_name', &
               'timestep count at the most recent accumulation, for each time-average stream')
       end if
    end if

    ! adding projection info
    if (glimmap_allocated(model%projection)) then
       status = parallel_def_var(NCO%id,glimmer_nc_mapvarname,NF90_CHAR,mapid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       call glimmap_CFPutProj(NCO%id,mapid,model%projection)
    end if

    ! setting the size of the level and staglevel dimension
    NCO%nlevel = model%general%upn
    NCO%nstaglevel = model%general%upn-1
    NCO%nstagwbndlevel = model%general%upn ! MJH this is the max index, not the size

    ! Note: The following dimension lengths are set to 1 by default,
    !       but can be increased depending on the config options.

    ! vertical coordinate for ocean data
    NCO%nzocn = model%ocean_data%nzocn

    ! vertical coordinate for atmosphere data
    NCO%nzatm = model%climate%nzatm

    ! glacier ID coordinate for glacier data
    NCO%nglacier = model%glacier%nglacier

    ! basin coordinate for basin data
    NCO%nbasin = model%ocean_data%nbasin

    ! coordinate for calvingMIP axis data
    NCO%naxis = model%calving%naxis

  end subroutine glimmer_nc_createfile

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_newfile(outfile, filename)

    !> Prepare an output object for writing to a new file, as for one_file_per_write.
    !> Call before glimmer_nc_createfile and NAME_io_create.
    !> Note: NAME_io_create consumes NCO%vars (and expands 'restart'), so the variable list
    !>       must be restored from NCO%vars_copy for each new file.
    !> Averaging state (processed_time, total_time, accum_tstep_count) is not changed,
    !>  so averages and time bounds carry over from one file to the next.

    implicit none
    type(glimmer_nc_output), pointer :: outfile
    character(len=*), intent(in) :: filename   !> name of the new file

    NCO%filename = filename
    NCO%vars = NCO%vars_copy
    outfile%timecounter = 1

  end subroutine glimmer_nc_newfile

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_closefile(outfile)

    !> Close the netCDF file for this output object, keeping the object itself.
    !> Used for one_file_per_write, where each write goes to a new file.

    use glimmer_log
    use glimmer_filenames
    use cism_parallel, only: parallel_close
    implicit none
    type(glimmer_nc_output), pointer :: outfile

    integer :: status

    if (NCO%file_open) then
       status = parallel_close(NCO%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
       NCO%file_open = .false.
       NCO%define_mode = .false.
       call write_log('Closing output file '//trim(process_path(NCO%filename)))
    end if

  end subroutine glimmer_nc_closefile

  !------------------------------------------------------------------------------

  function glimmer_nc_slice_filename(outfile, time) result(filename)

    !> Build the name of a single-slice file for one_file_per_write in standalone runs,
    !>  by inserting the time (yr) before the '.nc' suffix of outfile%base_filename.
    !> For example, 'out.tavg.nc' at time 6 becomes 'out.tavg.0006.nc'.
    !> Integer years are written with at least 4 digits (e.g., 0006, 1861, -21000);
    !>  other times are written with 3 decimal places (e.g., 0.500).
    !> Note: External drivers (e.g., the CESM wrapper) supply their own file names.

    implicit none
    type(glimmer_nc_output), pointer :: outfile
    real(dp), intent(in) :: time          ! model time (yr)
    character(len=fname_length) :: filename

    real(dp), parameter :: eps = 1.d-6    ! tolerance (yr) for treating time as an integer year
    character(len=fname_length) :: base
    character(len=32) :: time_str
    integer :: n

    base = outfile%base_filename
    if (len_trim(base) == 0) base = NCO%filename

    if (abs(time - real(nint(time),dp)) < eps) then
       write(time_str,'(i0.4)') nint(time)
    else
       write(time_str,'(f20.3)') time
       time_str = adjustl(time_str)
    end if

    n = len_trim(base)
    if (n > 3) then
       if (base(n-2:n) == '.nc') then
          filename = base(1:n-3) // '.' // trim(time_str) // '.nc'
          return
       end if
    end if
    filename = trim(base) // '.' // trim(time_str) // '.nc'

  end function glimmer_nc_slice_filename

  !------------------------------------------------------------------------------

  function glimmer_nc_output_has_var(outfile, varname) result(has_var)

    !> Return true if the variable varname belongs to this output object.
    !> If the file is open, check whether the file contains the variable.
    !> If no file is open (as for one_file_per_write, between writes), check the variable list.
    !> Note: The variable-list check does not expand the 'restart' keyword.
    !>       This is not a problem, since one_file_per_write is not allowed for restart files.

    implicit none
    type(glimmer_nc_output), pointer :: outfile
    character(len=*), intent(in) :: varname
    logical :: has_var

    integer :: status, varid

    if (NCO%file_open) then
       status = parallel_inq_varid(NCO%id, varname, varid)
       has_var = (status == NF90_NOERR)
    else
       has_var = (index(' '//trim(adjustl(NCO%vars_copy))//' ', ' '//trim(varname)//' ') /= 0)
    end if

  end function glimmer_nc_output_has_var

  !------------------------------------------------------------------------------

  function glimmer_nc_is_tavg_stream(outfile) result(is_tavg)

    !> Return true if this output object is a time-average stream, i.e., its variable list
    !>  (from the config file) includes time-average fields (names ending in '_tavg').
    !> The averaging state of these streams is saved in restart files.
    !> Note: Like the test for time bounds in glimmer_nc_createfile, this test uses the list from the
    !>       config file, so restart files (with the 'restart' keyword) are not time-average streams.

    implicit none
    type(glimmer_nc_output), pointer :: outfile
    logical :: is_tavg

    ! Note: Match whole words ending in '_tavg', so that '_tavg_sum' (restart files only) does not count.
    is_tavg = (index(trim(NCO%vars_copy)//' ', '_tavg ') /= 0)

  end function glimmer_nc_is_tavg_stream

  !------------------------------------------------------------------------------

  function glimmer_nc_is_restart(outfile) result(is_restart)

    !> Return true if this output object is a restart file, i.e., its variable list
    !>  (from the config file) includes the keyword 'restart'.

    implicit none
    type(glimmer_nc_output), pointer :: outfile
    logical :: is_restart

    is_restart = (index(' '//trim(adjustl(NCO%vars_copy))//' ', ' restart ') /= 0)

  end function glimmer_nc_is_restart

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_write_tavg_state(outfile, model)

    !> Write the averaging state of each time-average output stream (total time, start of the
    !>  averaging interval, timestep of the most recent accumulation) to a restart file,
    !>  at the current time slice.

    use glide_types
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    type(glide_global_type) :: model

    type(glimmer_nc_output), pointer :: p
    integer :: istream, status

    ! The streams are in the same order as in the global attribute 'tavg_streams' (see glimmer_nc_createfile).
    istream = 0
    p => model%funits%out_first
    do while (associated(p))
       if (glimmer_nc_is_tavg_stream(p)) then
          istream = istream + 1
          status = parallel_put_var(NCO%id, NCO%tavg_total_time_var, p%total_time, &
               (/istream, outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_put_var(NCO%id, NCO%tavg_start_time_var, p%nc%processed_time, &
               (/istream, outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_put_var(NCO%id, NCO%tavg_start_external_time_var, p%nc%processed_external_time, &
               (/istream, outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_put_var(NCO%id, NCO%tavg_accum_tstep_count_var, p%accum_tstep_count, &
               (/istream, outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
       end if
       p => p%next
    end do

  end subroutine glimmer_nc_write_tavg_state

  !------------------------------------------------------------------------------

  function glimmer_nc_stream_name(outfile) result(name)

    !> Return the name of an output stream, used to identify its averaging state in restart files:
    !>  base_filename (the name in the config file) for one_file_per_write streams
    !>  (e.g., 'h0a' in CESM), else the file name.

    implicit none
    type(glimmer_nc_output), pointer :: outfile
    character(len=fname_length) :: name

    if (len_trim(outfile%base_filename) > 0) then
       name = outfile%base_filename
    else
       name = NCO%filename
    end if

  end function glimmer_nc_stream_name

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_tavg_streams(model, nstreams, stream_list)

    !> Return the number of time-average output streams, and their names (space-separated, in the
    !>  order of the output list). This order is used for the averaging state in restart files.

    use glide_types
    implicit none
    type(glide_global_type) :: model
    integer, intent(out) :: nstreams
    character(len=*), intent(out) :: stream_list

    type(glimmer_nc_output), pointer :: p

    nstreams = 0
    stream_list = ''
    p => model%funits%out_first
    do while (associated(p))
       if (glimmer_nc_is_tavg_stream(p)) then
          nstreams = nstreams + 1
          stream_list = trim(stream_list)//' '//trim(glimmer_nc_stream_name(p))
       end if
       p => p%next
    end do
    stream_list = adjustl(stream_list)

  end subroutine glimmer_nc_tavg_streams

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_restore_tavg_state(outfile, model, restored)

    !> On a standard restart (restart = 1), restore the averaging state of a time-average output stream
    !>  (total_time, the start of the averaging interval in internal and external time, and
    !>  accum_tstep_count) from the restart file, so that the average continues exactly.
    !> The running sums (the '_tavg' arrays) are read separately, as '_tavg_sum' variables.
    !> The restart file is the first input file. If it holds no averaging state (e.g., it was written
    !>  by an older version of CISM), restored = .false. on return.
    !> The stream is found by name in the global attribute 'tavg_streams'. If the restart file has an
    !>  averaging state but not for this stream, the model aborts: the time-average output streams
    !>  should not change across a standard restart.

    use glimmer_log
    use glide_types
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    type(glide_global_type) :: model
    logical, intent(out) :: restored

    type(glimmer_nc_input), pointer :: infile
    character(len=tavg_list_len) :: stream_list, word
    character(len=fname_length) :: name
    character(len=msglen) :: message
    integer :: status, varid, istream, nstreams, i1, i2, n
    integer, dimension(2) :: start

    restored = .false.
    infile => model%funits%in_first
    if (.not. associated(infile)) return

    stream_list = ''
    status = parallel_get_att(NCI%id, NF90_GLOBAL, 'tavg_streams', stream_list)
    if (status /= NF90_NOERR) return   ! no averaging state in this file

    ! Find this stream in the list
    name = glimmer_nc_stream_name(outfile)
    istream = 0
    nstreams = 0
    n = len_trim(stream_list)
    i1 = 1
    do while (i1 <= n)
       if (stream_list(i1:i1) == ' ') then
          i1 = i1 + 1
          cycle
       end if
       i2 = i1 + index(stream_list(i1:)//' ', ' ') - 2   ! last character of this name
       nstreams = nstreams + 1
       word = stream_list(i1:i2)
       if (trim(word) == trim(name)) istream = nstreams
       i1 = i2 + 2
    end do

    if (istream == 0) then
       call write_log('Time-average output stream '//trim(name)//' has no averaging state in the restart file '// &
            trim(NCI%filename)//', which has averaging state for streams: '//trim(stream_list))
       call write_log('The time-average output streams should not change across a standard restart '// &
            '(restart = 1); to change them, use a hybrid restart (restart = 2)', GM_FATAL)
    end if

    ! Read the values for this stream at the time slice read from the restart file
    ! (dimensions: tavgstream, time)
    start = (/istream, infile%current_time/)

    status = parallel_inq_varid(NCI%id, 'tavg_total_time', varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_get_var(NCI%id, varid, outfile%total_time, start)
    call nc_errorhandle(__FILE__,__LINE__,status)

    status = parallel_inq_varid(NCI%id, 'tavg_start_time', varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_get_var(NCI%id, varid, NCO%processed_time, start)
    call nc_errorhandle(__FILE__,__LINE__,status)

    status = parallel_inq_varid(NCI%id, 'tavg_start_external_time', varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_get_var(NCI%id, varid, NCO%processed_external_time, start)
    call nc_errorhandle(__FILE__,__LINE__,status)

    status = parallel_inq_varid(NCI%id, 'tavg_accum_tstep_count', varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_get_var(NCI%id, varid, outfile%accum_tstep_count, start)
    call nc_errorhandle(__FILE__,__LINE__,status)

    restored = .true.

    write(message,*) 'Restored averaging state for time-average output ', trim(name), &
         ': interval start =', NCO%processed_time, ', total_time =', outfile%total_time
    call write_log(trim(message))

  end subroutine glimmer_nc_restore_tavg_state

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_checkwrite_init(outfile, time, external_time, tstep_count)

    !> Set the start of the first averaging interval for this output file.
    !> This is necessary to get the correct initial time bounds for tavg files
    !>  if this is a restart or a run not starting at t = 0.
    !> If tstep_count is present, then averages are not accumulated again until the model
    !>  takes a step beyond tstep_count.
    !> Called from openall_out when the file is created or reopened.
    !> Note: An external driver (e.g., the CESM wrapper) that passes an external time
    !>       to glimmer_nc_checkwrite should also call this subroutine with the external start time.

    implicit none
    type(glimmer_nc_output), pointer :: outfile
    real(dp), intent(in) :: time                     ! internal start time (yr)
    real(dp), intent(in), optional :: external_time  ! external start time; if not present, defaults to time
    integer, intent(in), optional :: tstep_count     ! timestep count at the start of the averaging interval

    NCO%processed_time = time
    if (present(external_time)) then
       NCO%processed_external_time = external_time
    else
       NCO%processed_external_time = time
    end if

    if (present(tstep_count)) then
       outfile%accum_tstep_count = tstep_count
    end if

  end subroutine glimmer_nc_checkwrite_init

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_checkwrite(outfile,model,forcewrite,time,external_time,wrote_timeslice)

    !> Check whether output is due for this file, and if so, write the time-slice variables.
    !> The caller (e.g., NAME_io_writeall) writes the model fields if wrote_timeslice = T.
    !>
    !> This subroutine is a driver for three steps, which can also be called separately
    !>  (e.g., by an external driver that controls when output is written):
    !> (1) glimmer_nc_advance_timeslice: leave define mode; advance the time counter after a write
    !> (2) glimmer_nc_write_due: decide whether a write is due
    !> (3) glimmer_nc_write_timeslice: write the time-slice variables

    use glimmer_log
    use glide_types
    use glimmer_filenames
    implicit none
    type(glimmer_nc_output), pointer :: outfile    
    type(glide_global_type) :: model
    logical forcewrite
    real(dp),optional :: time  ! time in years (written to 'internal_time')
    real(dp),optional :: external_time  ! external time (written to 'time'); if not present, defaults to internal_time
                                        ! units of external time are not necessarily years; e.g., CESM uses days
    logical, intent(out), optional :: wrote_timeslice  ! true if a time slice was written during this call
    real(dp) :: sub_time  ! local version of time (years)
    real(dp) :: sub_external_time  ! local version of external_time
    logical :: write_now  ! true if a time slice is written during this call

    ! Check for optional time argument
    if (present(time)) then
       sub_time = time
    else
       sub_time = model%numerics%time
    end if

    if (present(external_time)) then
       sub_external_time = external_time
    else
       sub_external_time = sub_time
    end if

    if (verbose_ncio .and. main_task) &
         write(iulog,*) 'In glimmer_nc_checkwrite, time, file =', sub_time, trim(process_path(NCO%filename))

    call glimmer_nc_advance_timeslice(outfile, sub_time)

    write_now = glimmer_nc_write_due(outfile, model, forcewrite, sub_time)

    if (write_now) then
       call glimmer_nc_write_timeslice(outfile, model, sub_time, sub_external_time)
    end if

    if (present(wrote_timeslice)) wrote_timeslice = write_now

  end subroutine glimmer_nc_checkwrite

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_advance_timeslice(outfile, time)

    !> Leave define mode if needed.
    !> If the file was written at an earlier time, then increment the time counter,
    !>  so that the next write goes to a new time slice.
    !> Note: This subroutine is needed for files with multiple time slices.

    implicit none
    type(glimmer_nc_output), pointer :: outfile
    real(dp), intent(in) :: time  ! current model time (yr)

    integer :: status

    ! Note: In glimmer_ncdf.F90, NCO%processed_time and NCO%processed_external_time
    !        are initialized to 0.0, and NCO%just_processed is initialized to FALSE.

    ! check if we are still in define mode and if so leave it
    ! Note: With one_file_per_write, there may be no open file.
    if (NCO%define_mode .and. NCO%file_open) then
       status = parallel_enddef(NCO%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
       NCO%define_mode = .FALSE.
    end if

    if (time > NCO%processed_time) then
       if (NCO%just_processed) then
          ! Finished writing during an earlier time step.
          ! For a file with multiple time slices, increase the counter so the next write goes to a new slice.
          ! With one_file_per_write, the file is already closed, and the next file will start at timecounter = 1.
          if (.not. outfile%one_file_per_write) then
             outfile%timecounter = outfile%timecounter + 1
             status = parallel_sync(NCO%id)
             call nc_errorhandle(__FILE__,__LINE__,status)
          end if
          NCO%just_processed = .FALSE.
       end if
    end if

  end subroutine glimmer_nc_advance_timeslice

  !------------------------------------------------------------------------------

  function glimmer_nc_write_due(outfile, model, forcewrite, time) result(write_due)

    !> Return true if output should be written to this file at the current time.
    !> This function sets no flags and writes nothing to the file, but it can write a warning to the log.

    use glimmer_log
    use glide_types
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    type(glide_global_type) :: model
    logical, intent(in) :: forcewrite  ! if true, write regardless of the output frequency
    real(dp), intent(in) :: time       ! current model time (yr)
    logical :: write_due

    character(len=msglen) :: message
    integer :: nfreq      ! freq/tinc; write output every nfreq timesteps
    real(dp) :: eps       ! tolerance for comparing time to end_write, to allow for roundoff

    write_due = .false.

    ! Compute the desired integer frequency for writing output (every nfreq timesteps), rounding if needed.
    ! Note: Both outfile%freq and model%general%tinc have units of years.
    !       If tinc does not divide evenly into freq, then output will be written at regular intervals,
    !        but not exactly at the user-desired frequency. For example, suppose tinc = 0.3 yr and freq = 1.0 yr.
    !       Then output will be written every 3 timesteps, since nint(1.0/0.3) = 3.
    !
    nfreq = nint(outfile%freq / model%numerics%tinc)

    if (nfreq == 0) then  ! freq < tinc/2
       nfreq = 1
       write(message,*) 'WARNING: output file frequency is smaller than timestep; writing output every timestep'
       call write_log(trim(message))
    endif

    ! Write output if any of the following is true:
    ! (1) forcewrite = T
    ! (2) tstep_count = 0 & write_init = T
    ! (3) tstep_count > 0 & mod(tstep_count,nfreq) = 0
    ! Note: write_init = T by default, but can be turned off in the config file (e.g., for restart files)
    ! In each case, the time must not be later than end_write, and the file must not already
    !  have been written at this time.
    ! Note: Model time is accumulated each timestep and can be slightly larger than end_write
    !       (e.g., 8.00000000000001 when end_write = 8), so allow a small tolerance.
    !       This follows the convention used for forcing times in NAME_read_forcing.

    eps = model%numerics%tinc * 1.0d-3

    if ( forcewrite .or.  &
        (model%numerics%tstep_count == 0 .and. outfile%write_init) .or.  &
        (model%numerics%tstep_count > 0 .and. mod(model%numerics%tstep_count, nfreq) == 0) ) then

       if (time <= outfile%end_write + eps .and. .not.NCO%just_processed) then
          write_due = .true.
       end if

    end if

  end function glimmer_nc_write_due

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_write_timeslice(outfile, model, time, external_time)

    !> Write internal_time, time, tstep_count, and (for tavg files) the time bounds
    !>  to the current time slice (outfile%timecounter).
    !> Then set NCO%just_processed = T, and set the start of the next averaging period.
    !>
    !> For instantaneous files, internal_time and time are the times passed in (the end of the
    !>  timestep). For time-average files (with time bounds), internal_time and time are the
    !>  midpoints of the averaging interval, following CESM and CF conventions; the bounds
    !>  give the start and end of the interval.

    use glimmer_log
    use glide_types
    use glimmer_filenames
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    type(glide_global_type) :: model
    real(dp), intent(in) :: time           ! time in years (written to 'internal_time')
    real(dp), intent(in) :: external_time  ! external time (written to 'time')

    character(len=msglen) :: message
    integer :: status
    real(dp) :: &
         internal_time_out,           & ! value written to internal_time
         external_time_out              ! value written to time
    real(dp), dimension(2) :: &
         internal_time_bounds,        & ! start and end times for averaging (internal)
         external_time_bounds           ! start and end times for averaging (external)

    ! Make sure the file is in data mode.
    ! Note: Normally this is done already in glimmer_nc_advance_timeslice.
    if (NCO%define_mode) then
       status = parallel_enddef(NCO%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
       NCO%define_mode = .FALSE.
    end if

    call write_log_div
    write(message,*) 'Writing to file ', trim(process_path(NCO%filename)), ' at time ', time
    call write_log(trim(message))

    if (verbose_ncio .and. main_task) &
         write(iulog,*) 'Writing to file ', trim(process_path(NCO%filename)), ' at time ', time

    ! Set the time values to write: the end of the timestep for instantaneous files,
    !  and the midpoint of the averaging interval for time-average files
    if (NCO%time_bounds) then
       internal_time_out = 0.5d0 * (NCO%processed_time + time)
       external_time_out = 0.5d0 * (NCO%processed_external_time + external_time)
    else
       internal_time_out = time
       external_time_out = external_time
    endif

    ! write time and tstep_count
    status = parallel_put_var(NCO%id, NCO%internal_timevar, internal_time_out, (/outfile%timecounter/))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_var(NCO%id, NCO%timevar, external_time_out, (/outfile%timecounter/))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_var(NCO%id, NCO%tstep_count_var, model%numerics%tstep_count, &
         (/outfile%timecounter/))
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! Write the time bounds if the file has them (i.e., a time-average file).
    ! Note: Test NCO%time_bounds rather than outfile%do_averages. A restart file can have
    !       do_averages = T (if its expanded variable list includes names with '_tavg'),
    !       but it has no time bounds variables.
    if (NCO%time_bounds) then

       write(message,*) '  Averaging interval (yr):', NCO%processed_time, time, ', total_time =', outfile%total_time
       call write_log(trim(message))

       internal_time_bounds(1) = NCO%processed_time
       internal_time_bounds(2) = time
       external_time_bounds(1) = NCO%processed_external_time
       external_time_bounds(2) = external_time

       if (verbose_ncio .and. main_task) &
            write(iulog,*) 'put internal_time_bounds:', internal_time_bounds(:)
       status = parallel_put_var(NCO%id, NCO%internal_timebounds_var, internal_time_bounds, &
            (/1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)

       if (verbose_ncio .and. main_task) &
            write(iulog,*) 'put external_time_bounds:', external_time_bounds(:)
       status = parallel_put_var(NCO%id, NCO%timebounds_var, external_time_bounds, &
            (/1,outfile%timecounter/))
       call nc_errorhandle(__FILE__,__LINE__,status)

    endif   ! time_bounds

    ! For restart files, write the averaging state of each time-average output stream
    if (NCO%tavg_state) call glimmer_nc_write_tavg_state(outfile, model)

    NCO%just_processed = .TRUE.

    ! reset the processed time for the next averaging period
    NCO%processed_time = time
    NCO%processed_external_time = external_time

  end subroutine glimmer_nc_write_timeslice

  !*****************************************************************************
  ! netCDF input
  !*****************************************************************************  

  subroutine openall_in(model)

    !> open all netCDF files for input
    use glide_types
    use glimmer_ncdf
    implicit none

    type(glide_global_type) :: model
    
    ! local variables
    type(glimmer_nc_input), pointer :: ic

    ! open input files
    ic=>model%funits%in_first
    do while(associated(ic))
       call glimmer_nc_openfile(ic,model)
       ic=>ic%next
    end do

    ! open forcing files
    ic=>model%funits%frc_first
    do while(associated(ic))
       call glimmer_nc_openfile(ic,model)
       ic=>ic%next
    end do

  end subroutine openall_in

  !------------------------------------------------------------------------------

  subroutine closeall_in(model)

    !> close all netCDF files for input
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    
    ! local variables
    type(glimmer_nc_input), pointer :: ic

    ! Input files
    ic=>model%funits%in_first
    do while(associated(ic))
       ic=>delete(ic)
    end do
    model%funits%in_first=>NULL()

    ! Forcing files
    ic=>model%funits%frc_first
    do while(associated(ic))
       ic=>delete(ic)
    end do
    model%funits%frc_first=>NULL()

  end subroutine closeall_in

  !------------------------------------------------------------------------------
  !TODO - Modify so the input file does not have to contain (x1,y1); OK if it just has (x0,y0)

  subroutine glimmer_nc_openfile(infile, model)

    !> open an existing netCDF file
    use glide_types
    use glimmer_map_CFproj
    use glimmer_map_types
    use glimmer_log
    use glimmer_filenames

    implicit none

    type(glimmer_nc_input), pointer :: infile    !> structure containg input netCDF descriptor
    type(glide_global_type) :: model             !> the model instance

    ! local variables
    integer dimsize, dimid, varid
    real, dimension(2) :: delta
    integer status    
    character(len=msglen) message
    
    real,parameter :: small = 1.e-6

    ! open netCDF file
    status = parallel_open(process_path(NCI%filename),NF90_NOWRITE,NCI%id)
    if (status /= NF90_NOERR) then
       call write_log('Error opening file '//trim(process_path(NCI%filename))//': '//nf90_strerror(status),&
            type=GM_FATAL,file=__FILE__,line=__LINE__)
    end if
    call write_log_div
    call write_log('opening file '//trim(process_path(NCI%filename))//' for input')

    ! getting projection, if none defined already
    if (.not.glimmap_allocated(model%projection)) model%projection = glimmap_CFGetProj(NCI%id)

    ! getting time dimension
    status = parallel_inq_dimid(NCI%id, 'time', NCI%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    ! get id of time variable
    status = parallel_inq_varid(NCI%id,glimmer_nc_internal_time_varname,NCI%internal_timevar)

    ! BACKWARDS_COMPATIBILITY(wjs, 2017-04-28) Older files may not have 'internal_time',
    ! so if we can't find that variable, fall back on 'time'.
    if (status /= NF90_NOERR) then
       status = parallel_inq_varid(NCI%id,glimmer_nc_time_varname,NCI%internal_timevar)
    end if
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! getting length of time dimension and allocating memory for array containing times
    status = parallel_inquire_dimension(NCI%id,NCI%timedim,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    allocate(infile%times(dimsize))
    infile%nt=dimsize
    status = parallel_get_var(NCI%id,NCI%internal_timevar,infile%times)

    ! getting tstep_count
    status = parallel_inq_varid(NCI%id,glimmer_nc_tstep_count_varname,NCI%tstep_count_var)
    ! BACKWARDS_COMPATIBILITY(wjs, 2017-05-17) Older files may not have 'tstep_count', so
    ! only read it if present.
    if (status == NF90_NOERR) then
       allocate(infile%tstep_counts(infile%nt))
       status = parallel_get_var(NCI%id,NCI%tstep_count_var,infile%tstep_counts)
       call nc_errorhandle(__FILE__,__LINE__,status)
       infile%tstep_counts_read = .true.
    else
       infile%tstep_counts_read = .false.
    end if

    ! setting the size of the level and staglevel dimension
    NCI%nlevel = model%general%upn
    NCI%nstaglevel = model%general%upn-1
    NCI%nstagwbndlevel = model%general%upn !MJH This is the max index, not size

    ! Note: The following dimension lengths are set to 1 by default,
    !       but can be increased depending on the config options.

    ! vertical coordinate for ocean data
    NCI%nzocn = model%ocean_data%nzocn

    ! vertical coordinate for atmosphere data
    NCI%nzatm = model%climate%nzatm

    ! glacier ID coordinate for glacier data
    NCI%nglacier = model%glacier%nglacier

    ! basin coordinate for basin data
    NCI%nbasin = model%ocean_data%nbasin

    ! coordinate for calvingMIP axis data
    NCI%naxis = model%calving%naxis

    ! checking if dimensions and grid spacing are the same as in the configuration file
    ! x1
    status = parallel_inq_dimid(NCI%id,'x1',dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize /= model%parallel%global_ewn) then
       write(message,*) 'Dimension x1 of file '//trim(process_path(NCI%filename))// &
            ' does not match with config dimension: ', dimsize, model%parallel%global_ewn
       call write_log(message,type=GM_FATAL)
    end if
    status = parallel_inq_varid(NCI%id,'x1',varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_get_var(NCI%id,varid,delta)
    call nc_errorhandle(__FILE__,__LINE__,status)

!WHL - mod to prevent code from crashing due to small roundoff error
!    if (abs(delta(2)-delta(1) - model%numerics%dew*len0) > small) then
    if (abs( (delta(2)-delta(1) - model%numerics%dew) / (model%numerics%dew) ) > small) then
       write(message,*) 'deltax1 of file '//trim(process_path(NCI%filename))// &
!!            ' does not match with config deltax: ', delta(2)-delta(1),model%numerics%dew*len0
            ' does not match with config deltax: ', delta(2)-delta(1),model%numerics%dew
       call write_log(message,type=GM_FATAL)
    end if

    ! x0
    !status = nf90_inq_dimid(NCI%id,'x0',dimid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_inquire_dimension(NCI%id,dimid,len=dimsize)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (dimsize /= model%general%ewn-1) then
    !   write(message,*) 'Dimension x0 of file ',trim(process_path(NCI%filename)),' does not match with config dimension: ', &
    !        dimsize, model%general%ewn-1
    !   call write_log(message,type=GM_FATAL)
    !end if
    !status = nf90_inq_varid(NCI%id,'x0',varid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_get_var(NCI%id,varid,delta)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (abs(delta(2)-delta(1) - model%numerics%dew) > small) then
    !   write(message,*) 'deltax0 of file '//trim(process_path(NCI%filename))//' does not match with config deltax: ', &
    !        delta(2)-delta(1),model%numerics%dew
    !   call write_log(message,type=GM_FATAL)
    !end if

    ! y1
    status = parallel_inq_dimid(NCI%id,'y1',dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize /= model%parallel%global_nsn) then
       write(message,*) 'Dimension y1 of file '//trim(process_path(NCI%filename))// &
            ' does not match with config dimension: ', dimsize, model%parallel%global_nsn
       call write_log(message,type=GM_FATAL)
    end if
    status = parallel_inq_varid(NCI%id,'y1',varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_get_var(NCI%id,varid,delta)
    call nc_errorhandle(__FILE__,__LINE__,status)


!WHL - mod to prevent code from crashing due to small roundoff error
!    if (abs(delta(2)-delta(1) - model%numerics%dns*len0) > small) then
    if (abs( (delta(2)-delta(1) - model%numerics%dns) / (model%numerics%dns) ) > small) then
       write(message,*) 'deltay1 of file '//trim(process_path(NCI%filename))// &
!!            ' does not match with config deltay: ', delta(2)-delta(1),model%numerics%dns*len0
            ' does not match with config deltay: ', delta(2)-delta(1),model%numerics%dns
       call write_log(message,type=GM_FATAL)
    end if
    
    ! y0
    !status = nf90_inq_dimid(NCI%id,'y0',dimid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_inquire_dimension(NCI%id,dimid,len=dimsize)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (dimsize /= model%general%nsn-1) then
    !   write(message,*) 'Dimension y0 of file '//trim(process_path(NCI%filename))//' does not match with config dimension: ',&
    !        dimsize, model%general%nsn-1
    !   call write_log(message,type=GM_FATAL)
    !end if
    !status = nf90_inq_varid(NCI%id,'y0',varid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_get_var(NCI%id,varid,delta)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (abs(delta(2)-delta(1) - model%numerics%dns) > small) then
    !   write(message,*) 'deltay0 of file '//trim(process_path(NCI%filename))//' does not match with config deltay: ',&
    !        delta(2)-delta(1),model%numerics%dns
    !   call write_log(message,type=GM_FATAL)
    !end if
  
  ! Check that the number of vertical layers is the same, though it's asking for trouble
  ! to check whether the spacing is the same (don't want to put that burden on setup,
  ! plus f.p. compare has been known to cause problems here)
  status = parallel_inq_dimid(NCI%id,'level',dimid)
  ! If we couldn't find the 'level' dimension, write a warning.
  ! We don't want to throw an error, as input files are only required to have it if they
  ! include 3D data fields.
  if (status == NF90_NOERR) then
        status = parallel_inquire_dimension(NCI%id, dimid, len=dimsize)
        call nc_errorhandle(__FILE__, __LINE__, status)
        if (dimsize /= model%general%upn .and. dimsize  /=  1) then
            write(message,*) 'Dimension level of file '//trim(process_path(NCI%filename))//&
                ' does not match with config dimension: ', &
                dimsize, model%general%upn
            call write_log(message,type=GM_FATAL)
        end if
  else
        call write_log("Input file contained no level dimension.  This is not necessarily a problem.", type=GM_WARNING)
  end if
  
  end subroutine glimmer_nc_openfile

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_checkread(infile,model,time)

    !> check if we should read from file
    !>
    !> If we're reading a restart file, then also sets model%numerics%tstart,
    !  model%numerics%time and model%numerics%tstep_count
    use glimmer_log
    use glide_types
    use glimmer_filenames

    implicit none

    type(glimmer_nc_input), pointer :: infile  !> structure containing output netCDF descriptor
    type(glide_global_type) :: model    !> the model instance
    real(dp),optional :: time           !> Optional alternative time

    character(len=msglen) :: message

    integer :: pos       ! to identify CISM standalone restart files
    integer :: pos_cesm  ! to identify CESM restart files

    real(dp) :: restart_time   ! time of restart (yr)

    ! Note: infile%nt = number of time slices in the file
    !       infile%current_time = current time index

    ! Parse the filename to see if it is a restart file (standalone or CESM)
    pos = index(infile%nc%filename,'.restart.') ! CISM naming convention for restart files
    pos_cesm = index(infile%nc%filename,'.r.')  ! CESM naming convention for restart files

    ! If a standalone file, then set current_time to the latest time slice
    ! (Not necessary for CESM restart files, which contain just one time slice)
    if (pos /= 0) then
       infile%current_time = infile%nt
    endif

    if (infile%current_time <= infile%nt) then
       if (.not.NCI%just_processed) then

          call write_log_div

          ! Reset model%numerics%tstart if reading a restart file
          !write(message,*) 'Check for restart:', trim(infile%nc%filename)
          !call write_log(message)

          if (pos /= 0 .or. pos_cesm /= 0) then   ! get the start time based on the current time slice

             if (model%options%is_restart == STANDARD_RESTART) then

                restart_time = infile%times(infile%current_time)      ! years
                model%numerics%tstart = restart_time
                model%numerics%time = restart_time

                if (infile%tstep_counts_read) then
                   model%numerics%tstep_count = infile%tstep_counts(infile%current_time)
                else
                   ! BACKWARDS_COMPATIBILITY(wjs, 2017-05-17) Older files may not have
                   ! 'tstep_count', so compute it ourselves here. We don't want to use this
                   ! formulation in general because it is prone to roundoff errors.
                   model%numerics%tstep_count = nint(model%numerics%time/model%numerics%tinc)
                end if

                write(message,*) 'Standard restart: New tstart, tstep_count =', &
                     model%numerics%tstart, model%numerics%tstep_count
                call write_log(message)

             elseif (model%options%is_restart == HYBRID_RESTART) then

                ! Use tstart from the config file, not the time from the restart file
                model%numerics%time = model%numerics%tstart  ! years
                model%numerics%tstep_count = 0

                write(message,*) 'Hybrid restart: New tstart, tstep_count =', &
                     model%numerics%tstart, model%numerics%tstep_count
                call write_log(message)

             endif   ! is_restart

          endif  ! pos/=0 or pos_cesm/=0

          write(message,*) 'Reading time slice ',infile%current_time,'(',infile%times(infile%current_time),') from file ', &
               trim(process_path(NCI%filename)), ' at time ', sub_time(model, time)
          call write_log(message)
          NCI%just_processed = .TRUE.
          NCI%processed_time = sub_time(model, time)

       end if  ! not just processed
    end if  ! current_time <= nt

    if (sub_time(model, time) > NCI%processed_time) then
       if (NCI%just_processed) then
          ! finished reading during last time step, need to increase counter...
          infile%current_time = infile%current_time + 1
          NCI%just_processed = .FALSE.
       end if
    end if

    ! For read_once files, suppress the call to glide_io_read by setting just_processed = false
    if (infile%read_once) then
       NCI%just_processed = .FALSE.
    endif

  contains

    real(dp) function sub_time(model, time)
      ! Get the current time applicable to this subroutine. 
      ! If time is present, use that; otherwise use model%numerics%time
      !
      ! We need this function to avoid code duplication. We canNOT simply set a local
      ! sub_time variable variable at the start of glimmer_nc_checkread, because model
      ! %numerics%time can be updated in the midst of this routine... so we need to
      ! determine sub_time when it's actually needed, with this function.
      use glide_types
      implicit none
      type(glide_global_type) :: model    !> the model instance
      real(dp),optional :: time           !> Optional alternative time

      if (present(time)) then
         sub_time = time
      else
         sub_time = model%numerics%time
      end if

    end function sub_time

  end subroutine glimmer_nc_checkread

  !------------------------------------------------------------------------------

  subroutine check_for_tempstag(whichdycore, nc)

      ! Check for the need to output tempstag and update the output variables if needed.
      !
      ! For the glam/glissade dycore, the vertical temperature grid has an extra level.
      ! In that case, the netCDF output file should include a variable
      ! called tempstag(0:nz) instead of temp(1:nz). This subroutine is added for
      ! convenience to allow the variable "temp" to be specified in the config
      ! file in all cases and have it converted to "tempstag" when appropriate.
      ! MJH
      !
      ! The same substitutions are applied to nc%vars and to nc%vars_copy.
      ! Note: Previously, this subroutine ended with nc%vars_copy = nc%vars. But NAME_io_create
      !       removes each variable from nc%vars as it is created, and this subroutine is called
      !       from NAME_io_create for each I/O module that writes to the file (e.g., glide_io and
      !       glad_io in CESM history files). For the second module, nc%vars had already been
      !       mostly consumed, so vars_copy lost most of the variable list. This matters when
      !       an output object writes more than one file (one_file_per_write), since each new file
      !       starts from vars_copy.

      use glimmer_log
      use glide_types

      implicit none
      integer, intent(in) :: whichdycore
      type(glimmer_nc_stat) :: nc

      ! If both temp and tempstag are specified, temp will get converted to tempstag
      ! and then there will be two tempstags in the list, but that is ok because
      ! the parser ignores duplicate entries in the varlist.
      ! (The check for the existence of variables looks like:    pos = index(NCO%vars,' acab ')  )

      ! Make sure vars_copy has a space at the beginning and end, as nc%vars does,
      ! so that the first and last variable names can be matched
      nc%vars_copy = ' '//trim(adjustl(nc%vars_copy))//' '

      if (whichdycore/=DYCORE_GLIDE) then
         ! We want temp, flwa and dissip to become tempstag, flwastag and dissipstag
         call replace_varname('temp', 'tempstag', &
              'Temperature remapping option uses temperature on a staggered vertical grid.' // &
              '  The netCDF output variable "temp" has been changed to "tempstag".')
         call replace_varname('flwa', 'flwastag', &
              'Temperature remapping option uses flwa on a staggered vertical grid.' // &
              '  The netCDF output variable "flwa" has been changed to "flwastag".')
         call replace_varname('dissip', 'dissipstag', &
              'Temperature remapping option uses dissip on a staggered vertical grid.' // &
              '  The netCDF output variable "dissip" has been changed to "dissipstag".')
      else  ! glide dycore
         ! We want tempstag, flwastag and dissipstag to become temp, flwa and dissip
         call replace_varname('tempstag', 'temp', &
              'The netCDF output variable "tempstag" should not be used with the Glide dycore.' // &
              '  The netCDF output variable "tempstag" has been changed to "temp".')
         call replace_varname('flwastag', 'flwa', &
              'The netCDF output variable "flwastag" should not be used with the Glide dycore.' // &
              '  The netCDF output variable "flwastag" has been changed to "flwa".')
         call replace_varname('dissipstag', 'dissip', &
              'The netCDF output variable "dissipstag" should not be used with the Glide dycore.' // &
              '  The netCDF output variable "dissipstag" has been changed to "dissip".')
      endif  ! whichdycore

    contains

      subroutine replace_varname(oldname, newname, message)

        ! Replace the first instance of variable oldname with newname, in both nc%vars and
        ! nc%vars_copy. If oldname is listed more than once, only the first instance is changed.
        ! Write the message to the log if nc%vars was changed.

        character(len=*), intent(in) :: oldname, newname, message
        logical :: changed

        call replace_word(nc%vars, oldname, newname, changed)
        if (changed) call write_log(message)
        call replace_word(nc%vars_copy, oldname, newname, changed)

      end subroutine replace_varname

      subroutine replace_word(str, oldname, newname, changed)

        ! Replace the first instance of ' oldname ' in str with ' newname '

        character(len=*), intent(inout) :: str
        character(len=*), intent(in) :: oldname, newname
        logical, intent(out) :: changed
        integer :: i

        i = index(str, ' '//oldname//' ')
        changed = (i > 0)
        if (changed) then
           str = str(1:i-1) // ' '//newname//' ' // str(i+len(oldname)+2:len(str))
        endif

      end subroutine replace_word

  end subroutine check_for_tempstag

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_get_dimlength(infile, dimname, dimlength)

    !WHL, Feb. 2022:
    ! This is a custom subroutine that opens an input file, finds the length
    ! of a specific dimension, and closes the file.
    ! It is useful for getting an array dimension whose size is not known in advance.
    ! Currently, it is called from glissade_initialise to get the length of the
    ! glacierid dimension, without having to put 'nglacier' in the config file by hand.

    use glimmer_ncdf
    use glimmer_log
    use glimmer_filenames, only: process_path

    type(glimmer_nc_input), pointer :: infile  !> structure containg input netCDF descriptor
    character(len=*), intent(in) :: dimname
    integer, intent(out) :: dimlength

    ! local variables
    integer :: status, dimid

    ! Open the file
    status = parallel_open(process_path(infile%nc%filename), NF90_NOWRITE, infile%nc%id)
    if (status /= NF90_NOERR) then
       call write_log('Error opening file '//trim(process_path(infile%nc%filename))//': '//nf90_strerror(status),&
            type=GM_FATAL, file=__FILE__,line=__LINE__)
    end if
    call write_log('Opening file '//trim(process_path(infile%nc%filename))//' for input')

    ! get the dimension length
    status = parallel_inq_dimid(infile%nc%id, trim(dimname), dimid)
    if (status .eq. nf90_noerr) then
       call write_log('Getting length of dimension '//trim(dimname)//' ')
       status = parallel_inquire_dimension(infile%nc%id, dimid, len=dimlength)
       if (status /= nf90_noerr) then
          call write_log('Error getting dimlength '//trim(dimname)//':'//nf90_strerror(status),&
               type=GM_FATAL, file=__FILE__,line=__LINE__)
       endif
    else
       call write_log('Error getting dimension '//trim(dimname)//':'//nf90_strerror(status),&
            type=GM_FATAL, file=__FILE__,line=__LINE__)
    endif

    ! close the file
    status = nf90_close(infile%nc%id)
    call write_log('Closing file '//trim(infile%nc%filename)//' ')

  end subroutine glimmer_nc_get_dimlength

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_get_var_integer_2d(infile, varname, field_2d)

    !WHL, July 2019:
    ! This is a custom subroutine that opens an input file, reads an integer array,
    ! and closes the file.  It is useful for reading fields that are needed for
    ! model initialization before the calls to openall_in and glide_io_readall.
    ! Currently, it is called from glissade_initialise to read ice_domain_mask,
    ! which is used to limit the computational domain to active blocks.

    use glimmer_ncdf
    use glimmer_log
    use glimmer_filenames, only: process_path

    type(glimmer_nc_input), pointer :: infile  !> structure containg input netCDF descriptor
    character(len=*), intent(in) :: varname
    integer, dimension(:,:), intent(inout) :: field_2d

    ! local variables
    integer :: status, varid

    ! Open the file
    status = parallel_open(process_path(infile%nc%filename), NF90_NOWRITE, infile%nc%id)
    if (status /= NF90_NOERR) then
       call write_log('Error opening file '//trim(process_path(infile%nc%filename))//': '//nf90_strerror(status),&
            type=GM_FATAL, file=__FILE__,line=__LINE__)
    end if
    call write_log('Opening file '//trim(process_path(infile%nc%filename))//' for input')

    ! read a 2D field
    status = parallel_inq_varid(infile%nc%id, trim(varname), varid)
    if (status .eq. nf90_noerr) then
       call write_log('Loading '//trim(varname)//' ')
       status = parallel_get_var(infile%nc%id, varid, field_2d)
       call nc_errorhandle(__FILE__,__LINE__, status)
    else
       call write_log('Error: Unable to read '//trim(varname)//' from file '//trim(process_path(infile%nc%filename))//' ', &
            type=GM_FATAL, file=__FILE__,line=__LINE__)
    endif

    ! close the file
    status = nf90_close(infile%nc%id)
    call write_log('Closing file '//trim(infile%nc%filename)//' ')

  end subroutine glimmer_nc_get_var_integer_2d

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_get_var_real8_2d(infile, varname, field_2d)

    !WHL, July 2019:
    ! This is the real8 version of the subroutine above.
    ! It is not currently called, but is included for generality.

    use glimmer_ncdf
    use glimmer_log
    use glimmer_filenames, only: process_path

    type(glimmer_nc_input), pointer :: infile  !> structure containg input netCDF descriptor
    character(len=*), intent(in) :: varname
    real(dp), dimension(:,:), intent(inout) :: field_2d

    ! local variables
    integer :: status, varid

    ! Open the file
    status = parallel_open(process_path(infile%nc%filename), NF90_NOWRITE, infile%nc%id)
    if (status /= NF90_NOERR) then
       call write_log('Error opening file '//trim(process_path(infile%nc%filename))//': '//nf90_strerror(status),&
            type=GM_FATAL, file=__FILE__,line=__LINE__)
    end if
    call write_log('Opening file '//trim(process_path(infile%nc%filename))//' for input')

    ! read a 2D field
    status = parallel_inq_varid(infile%nc%id, trim(varname), varid)
    if (status .eq. nf90_noerr) then
       call write_log('Loading '//trim(varname)//' ')
       status = parallel_get_var(infile%nc%id, varid, field_2d)
       call nc_errorhandle(__FILE__,__LINE__, status)
    else
       call write_log('Error: Unable to read '//trim(varname)//' from file '//trim(process_path(infile%nc%filename))//' ', &
            type=GM_FATAL, file=__FILE__,line=__LINE__)
    endif

    ! close the file
    status = nf90_close(infile%nc%id)
    call write_log('Closing file '//trim(infile%nc%filename)//' ')

  end subroutine glimmer_nc_get_var_real8_2d

!------------------------------------------------------------------------------


end module glimmer_ncio

!------------------------------------------------------------------------------
