!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_parallel_decomp.F90 - part of the Community Ice Sheet Model (CISM)  
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

module cism_parallel_decomp

  ! Part of the cism_parallel module:
  !  the domain decomposition (dividing the global grid among tasks) and the
  !  row and column communicators
  ! Other modules use these through cism_parallel.

  use glimmer_global, only : dp
  use glimmer_paramets, only: iulog
  use cism_parallel_core
  use cism_parallel_broadcast
  use cism_parallel_reduce
  use cism_parallel_gather_scatter

  implicit none

contains

!=======================================================================

  subroutine distributed_grid(ewn,      nsn,        &
                              parallel,             &
                              nhalo_in,             &
                              global_bc_in,         &
                              reprosum_in)

    ! Divide the global domain into blocks, with one task per block.
    ! Set various grid and domain variables for the local task.

    implicit none
    integer, intent(inout) :: ewn, nsn                  ! global grid dimensions
    type(parallel_type), intent(inout) :: parallel      ! info for parallel communication, computed here
    integer, intent(in), optional :: nhalo_in           ! number of rows of halo cells
    character(*), intent(in), optional :: global_bc_in  ! string indicating the global BC option
                                                        !TODO - Make this an integer?
    logical, intent(in), optional :: reprosum_in        ! if true, compute reproducible global sums

    integer :: best,i,j,metric
    real(dp) :: rewtasks,rnstasks
    integer :: ProcsEW

    ! begin

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         global_ewn  => parallel%global_ewn,   &
         global_nsn  => parallel%global_nsn,   &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn,      &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub,         &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         ewtasks     => parallel%ewtasks,      &
         nstasks     => parallel%nstasks,      &
         ewrank      => parallel%ewrank,       &
         nsrank      => parallel%nsrank,       &
         global_col_offset  => parallel%global_col_offset,  &
         global_row_offset  => parallel%global_row_offset,  &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub, &
         southwest_corner   => parallel%southwest_corner,   &
         southeast_corner   => parallel%southeast_corner,   &
         northeast_corner   => parallel%northeast_corner,   &
         northwest_corner   => parallel%northwest_corner,   &
         staggered_ilo      => parallel%staggered_ilo,      &
         staggered_ihi      => parallel%staggered_ihi,      &
         staggered_jlo      => parallel%staggered_jlo,      &
         staggered_jhi      => parallel%staggered_jhi,      &
         reprosum           => parallel%reprosum            &
         )

    ! set the boundary conditions (periodic by default)
    ! Note: The no-penetration BC is treated as periodic. This BC may need some more work.

    if (present(global_bc_in)) then
       if (trim(global_bc_in) == 'periodic') then
          periodic_bc = .true.
          outflow_bc = .false.
          no_ice_bc = .false.
          if (main_task) write(*,*) 'Setting periodic boundary conditions'
       elseif (trim(global_bc_in) == 'outflow') then
          periodic_bc = .false.
          outflow_bc = .true.
          no_ice_bc = .false.
          if (main_task) write(*,*) 'Setting outflow boundary conditions'
       elseif (trim(global_bc_in) == 'no_penetration') then
          periodic_bc = .true.   ! Currently use the same halo logic for no_penetration and periodic  
          outflow_bc = .false.
          no_ice_bc = .false.
          if (main_task) write(*,*) 'Setting no_penetration boundary conditions'
       elseif (trim(global_bc_in) == 'no_ice') then
          periodic_bc = .false.
          outflow_bc = .false.
          no_ice_bc = .true.
          if (main_task) write(*,*) 'Setting no_ice boundary conditions'
       else
          if (main_task) write(*,*) 'Error: Invalid global_bc option for distributed_grid subroutine'
          call parallel_stop(__FILE__, __LINE__)
       endif
    else   ! default to periodic
       periodic_bc = .true.
       outflow_bc = .false.
       no_ice_bc = .false.
    endif

    ! Optionally, change the halo values
    ! Note: The Glissade higher-order dycore requires nhalo = 2.
    !       The Glide SIA dycore requires nhalo = 0.
    ! The default halo values at the top of the module are appropriate for
    !  the higher-order dycores.  Here they can be reset to zero for Glide.

    if (present(nhalo_in)) then
       if (main_task) then
          write(*,*) 'Setting halo values: nhalo =', nhalo_in
          if (nhalo_in < 0) then
             write(*,*) 'ERROR: nhalo must be >= 0'
             call parallel_stop(__FILE__, __LINE__)
          elseif (nhalo_in /= 2) then
             write(*,*) 'WARNING: parallel dycores tested only with nhalo = 2'
          endif
       endif 
       nhalo = nhalo_in
       lhalo = nhalo
       uhalo = nhalo
       staggered_lhalo = lhalo
       staggered_uhalo = max(uhalo-1, 0)
    endif

    global_ewn = ewn
    global_nsn = nsn

    ewtasks = 0
    nstasks = 0
    best = huge(best)
    do i = 1,min(tasks,global_ewn)
       j = tasks/i
       if (j<=global_nsn.and.i*j==tasks) then ! try to use all tasks
          metric = abs(i*global_nsn-j*global_ewn) ! zero if ewn/nsn == i/j
          if (metric<best) then
             best = metric
             ewtasks = i
             nstasks = j
          end if
       end if
    end do
    if (ewtasks*nstasks/=tasks) call parallel_stop(__FILE__,__LINE__)

    ! Store critical value for creating global IDs.  Defines grid distribution.
    ProcsEW = ewtasks

    ! For globalID calculations, determine processor's global grid index offsets.
    ! Sum block sizes for row blocks preceding this_rank.
    ! Do not include halo offsets in global calculations.
    ! (There are ProcsEW processors per row.)
    global_col_offset = 0
    do ewrank=0,mod(this_rank, ProcsEW)-1
      rewtasks = 1/real(ewtasks,dp)
      ewlb = nint(ewrank*global_ewn*rewtasks)+1
      ewub = nint((ewrank+1)*global_ewn*rewtasks)
      own_ewn = ewub-ewlb+1
      global_col_offset = global_col_offset + own_ewn
    enddo

    ! Sum block sizes for column blocks preceding this_rank
    ! (Integer division required for this_rank/ProcsEW)
    global_row_offset = 0
    do nsrank=0,(this_rank/ProcsEW)-1
      rnstasks = 1/real(nstasks,dp)
      nslb = nint(nsrank*global_nsn*rnstasks)+1
      nsub = nint((nsrank+1)*global_nsn*rnstasks)
      own_nsn = nsub-nslb+1
      global_row_offset = global_row_offset + own_nsn
    enddo

    ! Set local processor's grid indices, including halo offsets
    ewrank = mod(this_rank,ewtasks)
    rewtasks = 1/real(ewtasks,dp)
    ewlb = nint(ewrank*global_ewn*rewtasks)+1-lhalo
    ewub = nint((ewrank+1)*global_ewn*rewtasks)+uhalo
    local_ewn = ewub-ewlb+1
    own_ewn = local_ewn-lhalo-uhalo
    ewn = local_ewn

    nsrank = this_rank/ewtasks
    rnstasks = 1/real(nstasks,dp)
    nslb = nint(nsrank*global_nsn*rnstasks)+1-lhalo
    nsub = nint((nsrank+1)*global_nsn*rnstasks)+uhalo
    local_nsn = nsub-nslb+1
    own_nsn = local_nsn-lhalo-uhalo
    nsn = local_nsn

    ! Determine the global bounds of ewlb, ewub, nslb, and nsub.
    ! These are used to allocate global arrays used in gathers and scatters.
    global_minval_ewlb = parallel_reduce_min(ewlb)
    global_maxval_ewub = parallel_reduce_max(ewub)
    global_minval_nslb = parallel_reduce_min(nslb)
    global_maxval_nsub = parallel_reduce_max(nsub)

    west = this_rank-1
    if ((west/ewtasks<this_rank/ewtasks).or.(west<0)) west = west+ewtasks
    east = this_rank+1
    if (east/ewtasks>this_rank/ewtasks) east = east-ewtasks
    south = this_rank-ewtasks
    if (south<0) south = south+tasks
    north = this_rank+ewtasks
    if (north>=tasks) north = north-tasks

    ! Specify that this is not a corner task.
    ! A southeast corner task is defined as a task with active south and east neighbors but an
    !  inactive southeast neighbor; and similarly for other directions.
    ! This subroutine assumes a standard rectangular layout of tasks, in which case
    !  there are no such corner tasks.
    ! Subroutine distributed_grid_active_blocks (below) typically distributes tasks
    !  such that there are corner tasks.
    southwest_corner = .false.
    southeast_corner = .false.
    northeast_corner = .false.
    northwest_corner = .false.

    ! Set the limits of locally owned vertices on the staggered grid.
    ! For periodic BC, staggered_ilo = staggered_jlo = staggered_lhalo+1.
    ! For outflow BC the locally owned vertices include the southern and western rows
    !  of the global domain, so staggered_ilo = staggered_jlo = staggered_lhalo on
    !  processors that include these rows.
    ! Note: For no_ice BC, we assume (uvel,vvel) = 0 along the global boundary.
    !       In this case, vertices along the southern and western edges of the global boundary
    !        are not considered to be locally owned by any task.

    if (outflow_bc .and. this_rank <= west) then  ! on west edge of global domain
       staggered_ilo = staggered_lhalo
    else
       staggered_ilo = staggered_lhalo+1
    endif
    staggered_ihi = ewn - 1 - staggered_uhalo

    if (outflow_bc .and. this_rank <= south) then  ! on south edge of global domain
       staggered_jlo = staggered_lhalo
    else
       staggered_jlo = staggered_lhalo+1
    endif
    staggered_jhi = nsn - 1 - staggered_uhalo
    
    ! Check that we have not split up the problem too much.  We do not want halos overlapping in either dimension.
    ! local_* - lhalo - uhalo is the actual number of non-halo cells on a processor.
    if ((local_nsn - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
        write(*,*) "NS halos overlap on processor ", this_rank
        call parallel_stop(__FILE__, __LINE__)
    endif

    if ((local_ewn  - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
        write(*,*) "EW halos overlap on processor ", this_rank
        call parallel_stop(__FILE__, __LINE__)
    endif

    if (present(reprosum_in)) then
       reprosum = reprosum_in
    else
       reprosum = .false.
    endif

    ! If computing reproducible sums, then set some options
    ! Note: For standalone CISM, reprosum = F by default; can set = T in the config file
    !       For CESM coupled runs, reprosum = T by default
    if (reprosum) then

    !WHL - Commented out the reprosum_setops calls.
    !      For standalone CISM, the call to cism_reprosum_setops is not needed.
    !      TBD whether this is true for coupled CESM runs.
#ifdef CCSM_COUPLED
!!       call shr_reprosum_setops()
#else
!!       call cism_reprosum_setops()
#endif
    endif  ! reprosum

!    call parallel_barrier
!    write(iulog,*) 'task, west, east, south, north:', this_rank, west, east, south, north

    ! Uncomment to print grid geometry
!    write(iulog,*) "Process ", this_rank, " Total = ", tasks, " ewtasks = ", ewtasks, " nstasks = ", nstasks
!    write(iulog,*) "Process ", this_rank, " ewrank = ", ewrank, " nsrank = ", nsrank
!    write(iulog,*) "Process ", this_rank, " l_ewn = ", local_ewn, " o_ewn = ", own_ewn
!    write(iulog,*) "Process ", this_rank, " l_nsn = ", local_nsn, " o_nsn = ", own_nsn
!    write(iulog,*) "Process ", this_rank, " ewlb = ", ewlb, " ewub = ", ewub
!    write(iulog,*) "Process ", this_rank, " nslb = ", nslb, " nsub = ", nsub
!    write(iulog,*) "Process ", this_rank, " east = ", east, " west = ", west
!    write(iulog,*) "Process ", this_rank, " north = ", north, " south = ", south
!    write(iulog,*) "Process ", this_rank, " ew_vars = ", own_ewn, " ns_vars = ", own_nsn
!    write(iulog,*) "Process ", this_rank, " global_col_offset = ", global_col_offset, &
!                                      " global_row_offset = ", global_row_offset

    call distributed_print_grid(own_ewn, own_nsn)

    end associate

  end subroutine distributed_grid

!=======================================================================

  subroutine distributed_grid_active_blocks(ewn,      nsn,        &
                                            nx_block, ny_block,   &
                                            ice_domain_mask,      &
                                            parallel,             &
                                            reprosum_in,          &
                                            inquire_only)

    ! Divide the global domain into blocks, setting various grid and domain variables
    !  for each block as in subroutine distributed_grid above.
    ! Then read a mask to identify which blocks are potentially active
    !  (i.e, ice can be present in some or all of the block).
    !  Assign a task to each active block.  Activate additional blocks as needed
    !  so that there is exactly one task per active block.
    ! Assume no_ice boundary conditions.  For these BCs, scalars (including ice thickness)
    !  are set to zero not only in the global halo (as for outflow BCs), but also
    !  along one row just inside the global boundary.  This ensures that velocity = 0
    !  for vertices along the global boundary, and allows halo routines to work
    !  correctly when some tasks are inactive.
    ! Set the neighbor task indices (west, east, south, north) such that if a task
    !  does not have an active west neighbor, it has west = this_rank (and similarly
    !  for other directions).
    ! Set logical variables for corner tasks.  A southeast corner task is defined as
    !  a task with active south and east neighbors but an inactive southeast neighbor
    !  (and similarly for other directions).
    !
    ! Note: The code does not yet compute the total ice mass zeroed out near the global boundary.
    !       This does not result in a mass conservation error, because the total dmass/dt term
    !        is computed by summing over dH/dt in locally owned cells only.
    !       There is a similar issue with outflow BC; the code does not compute the ice mass
    !        zeroed out in the global halo.
    ! TODO: Compute dmass/dt by summing over the global ice mass each time step?
    !       Then compute an outflow flux to bring mass conservation back into balance.
    !       This would require computing the total ice mass before and after halo updates.x

    ! Created by WHL, July 2019, based on subroutine distributed_grid above.

    implicit none

    integer, intent(inout) :: ewn, nsn              ! global grid dimensions
    integer, intent(in) :: nx_block, ny_block       ! block sizes in each direction
    integer, intent(in), dimension(:,:) :: &
         ice_domain_mask                            ! = 1 where ice is potentially present and active, else = 0
    type(parallel_type), intent(inout) :: parallel  ! info for parallel communication, computed here
    logical, intent(in), optional :: reprosum_in    ! if true, compute reproducible global sums
    logical, intent(in), optional :: inquire_only   ! if true, then report the number of active blocks and abort

    integer :: i, j, nb, nt
    integer :: nblocks               ! number of blocks = ewtasks * nstasks
    real(dp) :: rewtasks, rnstasks
    integer :: ProcsEW

    integer :: nblocks_active        ! number of active blocks
    logical :: only_inquire          ! local version of inquire_only
    integer :: corner_block          ! block number for a corner neighbor

    ! arrays with dimension 'nblocks'

    integer, dimension(:), allocatable ::  &
         ewrank_block,             &  ! grid column in which a block lies, increasing from W to E
         nsrank_block,             &  ! grid row in which a block lies, increasing from S to N
         ewlb_block, ewub_block,   &  ! lower and upper bounds in E-W direction for each block, including halos
         nslb_block, nsub_block,   &  ! lower and upper bounds in N-S direction for each block, including halos
         own_ewn_block,            &  ! number of cells in E-W direction on each block   !TODO - rename?  Not "owned"
         own_nsn_block,            &  ! number of cells in N-S direction on each block   !TODO - rename?  Not "owned"
         global_col_offset_block,  &  ! column offset for each block
         global_row_offset_block,  &  ! row offset for each block
         local_ewn_block,          &  ! number of cells in E-W direction on each block, including halos
         local_nsn_block,          &  ! number of cells in N-S direction on each block, including halos
         ewn_block,                &  ! set to local_ewn_block
         nsn_block,                &  ! set to local_nsn_block
         east_block, west_block,   &  ! east and west neighbors of each block
         north_block, south_block     ! north and south neighbors of each block

    integer, dimension(:), allocatable ::  &
         block_to_task                ! task (if any) assigned to each block

    logical, dimension(:), allocatable :: &
         block_is_active,          &  ! true for active blocks
         block_is_active_new          ! augmented version of block_is_active

    ! arrays with dimension 'tasks'

    integer, dimension(:), allocatable ::  &
         task_to_block                ! block associated with each task

    logical :: verbose_active_blocks = .false.

    associate(  &
         periodic_bc => parallel%periodic_bc,  &
         outflow_bc  => parallel%outflow_bc,   &
         no_ice_bc   => parallel%no_ice_bc,    &
         global_ewn  => parallel%global_ewn,   &
         global_nsn  => parallel%global_nsn,   &
         local_ewn   => parallel%local_ewn,    &
         local_nsn   => parallel%local_nsn,    &
         own_ewn     => parallel%own_ewn,      &
         own_nsn     => parallel%own_nsn,      &
         ewlb        => parallel%ewlb,         &
         ewub        => parallel%ewub,         &
         nslb        => parallel%nslb,         &
         nsub        => parallel%nsub,         &
         east        => parallel%east,         &
         west        => parallel%west,         &
         north       => parallel%north,        &
         south       => parallel%south,        &
         ewtasks     => parallel%ewtasks,      &
         nstasks     => parallel%nstasks,      &
         ewrank      => parallel%ewrank,       &
         nsrank      => parallel%nsrank,       &
         global_col_offset  => parallel%global_col_offset,  &
         global_row_offset  => parallel%global_row_offset,  &
         global_minval_ewlb => parallel%global_minval_ewlb, &
         global_maxval_ewub => parallel%global_maxval_ewub, &
         global_minval_nslb => parallel%global_minval_nslb, &
         global_maxval_nsub => parallel%global_maxval_nsub, &
         southwest_corner   => parallel%southwest_corner,   &
         southeast_corner   => parallel%southeast_corner,   &
         northeast_corner   => parallel%northeast_corner,   &
         northwest_corner   => parallel%northwest_corner,   &
         staggered_ilo      => parallel%staggered_ilo,      &
         staggered_ihi      => parallel%staggered_ihi,      &
         staggered_jlo      => parallel%staggered_jlo,      &
         staggered_jhi      => parallel%staggered_jhi,      &
         reprosum           => parallel%reprosum            &
         )

    if (present(inquire_only)) then
       only_inquire = inquire_only
       if (only_inquire) verbose_active_blocks = .true.
    else
       only_inquire = .false.
    endif

    ! Set the boundary conditions.
    ! For the active_blocks option, only the no_ice BC is supported.

    no_ice_bc = .true.
    outflow_bc = .false.
    periodic_bc = .false.

    if (main_task) write(iulog,*) 'Setting no_ice boundary conditions, assigning tasks to active blocks only'

    ! Note: There is no option to change nhalo.
    !       The Glissade higher-order dycore requires nhalo = 2, which is set with other halo values at the top of this module.
    !       The Glide SIA dycore (with nhalo = 0) is serial-only and is not compatible with the active-block option.

    global_ewn = ewn
    global_nsn = nsn

    ! Given the global dimensions global_ewn and global_nsn, along with the specified block size nx_block, ny_block,
    ! determine the number of blocks needed to span the global domain in each direction.
    ! Note: I considered renaming ewtasks and nstasks to ewblocks and nsblocks, but kept the old names
    !       because they're used elsewhere in the code.  With the active_block logic, the domain of active tasks
    !       might not fully span the region defined by (ewtasks,nstasks).

    if (nx_block == 0 .or. ny_block == 0) then
       write(iulog,*) 'Error: Must have nx_block and ny_block > 0'
       call parallel_stop(__FILE__, __LINE__)
    endif

    if (mod(global_ewn, nx_block) > 0) then
       ewtasks = global_ewn/nx_block + 1
    else
       ewtasks = global_ewn/nx_block
    endif

    if (mod(global_nsn, ny_block) > 0) then
       nstasks = global_nsn/ny_block + 1
    else
       nstasks = global_nsn/ny_block
    endif

    nblocks = ewtasks * nstasks

    if (main_task .and. verbose_active_blocks) then
       write(iulog,*) 'global_ewn, global_nsn =', global_ewn, global_nsn
       write(iulog,*) 'nx_block, ny_block =', nx_block, ny_block
       write(iulog,*) 'ewtasks, nstasks, nblocks =', ewtasks, nstasks, nblocks
    endif

    ! Allocate block variables
    allocate(ewrank_block(0:nblocks-1))
    allocate(global_col_offset_block(0:nblocks-1))
    allocate(ewlb_block(0:nblocks-1))
    allocate(ewub_block(0:nblocks-1))
    allocate(local_ewn_block(0:nblocks-1))
    allocate(own_ewn_block(0:nblocks-1))
    allocate(ewn_block(0:nblocks-1))

    allocate(nsrank_block(0:nblocks-1))
    allocate(global_row_offset_block(0:nblocks-1))
    allocate(nslb_block(0:nblocks-1))
    allocate(nsub_block(0:nblocks-1))
    allocate(local_nsn_block(0:nblocks-1))
    allocate(own_nsn_block(0:nblocks-1))
    allocate(nsn_block(0:nblocks-1))

    allocate(block_is_active(0:nblocks-1))
    allocate(block_is_active_new(0:nblocks-1))

    allocate(west_block(0:nblocks-1))
    allocate(east_block(0:nblocks-1))
    allocate(south_block(0:nblocks-1))
    allocate(north_block(0:nblocks-1))

    ! Determine properties of each block.
    ! Since not all blocks are active, we may not have a one-to-one correspondence between
    !  blocks and tasks, so we cannot yet compute task-specific properties on the local task.
    ! Note: These calculations could potentially be done on main_task only, then broadcast to local tasks.
    !       Doing them on all tasks avoids the need to broadcast arrays of size (nblocks) to all local tasks.

    ! For globalID calculations, determine each block's global grid index offsets.
    ! Do not include halo offsets in global calculations.

    ! Store critical value for creating global IDs.  Defines grid distribution.
    !Note: Currently, ProcsEW is used only in cism_profile.F90; could simply replace with ewtasks?
    ProcsEW = ewtasks

    ! Loop over blocks.  Note zero-based indexing.
    do nb = 0, nblocks-1

       ! Sum block sizes for row blocks preceding each block.
       ! The number of blocks per row is ewtasks.

       global_col_offset_block(nb) = 0
       do ewrank = 0, mod(nb,ewtasks) - 1
          rewtasks = 1.0d0/real(ewtasks,dp)
          ewlb_block(nb) = nint( ewrank    * global_ewn * rewtasks) + 1
          ewub_block(nb) = nint((ewrank+1) * global_ewn * rewtasks)
          own_ewn_block(nb) = ewub_block(nb) - ewlb_block(nb) + 1
          global_col_offset_block(nb) = global_col_offset_block(nb) + own_ewn_block(nb)
       enddo

       ! Sum block sizes for column blocks preceding each block.
       ! The number of blocks per column is nstasks.
       ! (Integer division required for nb/ewtasks)

       global_row_offset_block(nb) = 0
       do nsrank = 0, (nb/ewtasks) - 1
          rnstasks = 1.0d0/real(nstasks,dp)
          nslb_block(nb) = nint( nsrank * global_nsn * rnstasks) + 1
          nsub_block(nb) = nint((nsrank+1) * global_nsn * rnstasks)
          own_nsn_block(nb) = nsub_block(nb) - nslb_block(nb) + 1
          global_row_offset_block(nb) = global_row_offset_block(nb) + own_nsn_block(nb)
       enddo

       ! Set each block's grid indices, including halo offsets

       ewrank_block(nb) = mod(nb, ewtasks)
       rewtasks = 1.0d0/real(ewtasks,dp)
       ewlb_block(nb) = nint( ewrank * global_ewn * rewtasks) + 1 - lhalo
       ewub_block(nb) = nint((ewrank+1) * global_ewn * rewtasks) + uhalo
       local_ewn_block(nb) = ewub_block(nb) - ewlb_block(nb) + 1
       own_ewn_block(nb) = local_ewn_block(nb) - lhalo - uhalo
       ewn_block(nb) = local_ewn_block(nb)

       nsrank_block(nb) = nb / ewtasks
       rnstasks = 1.0d0/real(nstasks,dp)
       nslb_block(nb) = nint( nsrank * global_nsn * rnstasks) + 1 - lhalo
       nsub_block(nb) = nint((nsrank+1) * global_nsn * rnstasks) + uhalo
       local_nsn_block(nb) = nsub_block(nb) - nslb_block(nb) + 1
       own_nsn_block(nb) = local_nsn_block(nb) - lhalo - uhalo
       nsn_block(nb) = local_nsn_block(nb)

       ! Identify the west, east, south, and north neighbors of each block
       !
       ! Periodic BCs satisfy the convention that for blocks on the western edge of the domain,
       !  'west' is the easternmost block of that row, and for blocks on the easter edge of the domain,
       ! 'east' is the westernmost block of that row.  And similiarly for south and north.
       !
       ! For no_ice BC, we have a different convention: For blocks along the western edge of the domain
       !  (where west_block/ewtasks < nb/ewtasks by integer division), west_block is set to nb.
       ! According to the halo logic, when this_rank <= west (and specifically when this_rank = west),
       !  no message is sent to the west, and none will be received from the west.
       ! And similarly for other directions.

       west_block(nb) = nb - 1
       if ((west_block(nb)/ewtasks < nb/ewtasks) .or. (west_block(nb) < 0)) west_block(nb) = nb

       east_block(nb) = nb + 1
       if (east_block(nb)/ewtasks > nb/ewtasks) east_block(nb) = nb

       south_block(nb) = nb - ewtasks
       if (south_block(nb) < 0) south_block(nb) = nb

       north_block(nb) = nb + ewtasks
       if (north_block(nb) >= nblocks) north_block(nb) = nb

    enddo   ! nblocks

    ! Determine the global bounds of ewlb, ewub, nslb, and nsub.
    ! These are used to allocate global arrays used in gathers and scatters.
    global_minval_ewlb = minval(ewlb_block)
    global_maxval_ewub = maxval(ewub_block)
    global_minval_nslb = minval(nslb_block)
    global_maxval_nsub = maxval(nsub_block)

    !WHL - debug
    if (main_task) then
       write(iulog,*) 'global_minval_ewlb =', global_minval_ewlb
       write(iulog,*) 'global_maxval_ewub =', global_maxval_ewub
       write(iulog,*) 'global_minval_nslb =', global_minval_nslb
       write(iulog,*) 'global_maxval_nsub =', global_maxval_nsub
    endif

    ! Determine which blocks are active.
    ! A block is active if one or more locally owned cells has ice_domain_mask = 1.
    ! Note: This calculation is done on main_task only, since this is where
    !        the global array ice_domain_mask is nonzero.
    !       The following loop could be expensive on fine global grids.

    block_is_active(:) = .false.
    nblocks_active = 0

    if (main_task) then
       do nb = 0, nblocks-1
          ij_outer: do j = nslb_block(nb) + lhalo, nsub_block(nb) - uhalo
             do i = ewlb_block(nb) + lhalo, ewub_block(nb) - uhalo
                if (ice_domain_mask(i,j) == 1) then
                   block_is_active(nb) = .true.
                   nblocks_active = nblocks_active + 1
                   exit ij_outer
                endif
             enddo
          enddo ij_outer
       enddo

       if (verbose_active_blocks) then
          write(iulog,*) 'nblocks, nblocks_active:', nblocks, nblocks_active
          write(iulog,*) ' '
          write(iulog,*) 'Block layout:'
          do j = nstasks-1, 0, -1
             do i = 0, ewtasks-1
                nb = ewtasks*j + i
                write(iulog, '(i5)', advance='no') nb
             enddo
             write(iulog,*) ' '
          enddo
          write(iulog,*) ' '
          write(iulog,*) 'block_is_active:'
          do j = nstasks-1, 0, -1
             do i = 0, ewtasks-1
                nb = ewtasks*j + i
                write(iulog, '(l5)', advance='no') block_is_active(nb)
             enddo
             write(iulog,*) ' '
          enddo
       endif  ! verbose_active_blocks

       if (only_inquire) then  ! report the number of active blocks, then abort cleanly
          write(iulog,*)
          write(iulog,*) 'The number of active blocks with this domain and block layout is ', nblocks_active
          write(iulog,*) 'The total number of blocks is ', nblocks
          write(iulog,*) 'Please resubmit with nblocks_active <= tasks <= n_blocks'
          call parallel_stop(__FILE__, __LINE__)
       else   ! abort if tasks < nblocks_active or tasks > nblocks; otherwise proceed
          if (tasks < nblocks_active) then
             write(iulog,*)
             write(iulog,*) 'Fatal error: tasks < nblocks_active'
             write(iulog,*) 'Number of tasks =', tasks
             write(iulog,*) 'Minimum number of tasks to compute on all active blocks is ', nblocks_active
             call parallel_stop(__FILE__, __LINE__)
          elseif (tasks > nblocks) then
             write(iulog,*)
             write(iulog,*) 'Fatal error: tasks > nblocks'
             write(iulog,*) 'Number of tasks =', tasks
             write(iulog,*) 'Maximum number of tasks to compute on all blocks is ', nblocks
             call parallel_stop(__FILE__, __LINE__)
          endif
       endif   ! only_inquire

       ! If nblocks_active < tasks, then activate additional blocks along the ice sheet boundary until
       !  nblocks_active = tasks. Do not activate "orphan" blocks disconnected from other active blocks.

       if (nblocks_active < tasks) then

          if (verbose_active_blocks) then
             write(iulog,*) 'nblocks_active, tasks =', nblocks_active, tasks
             write(iulog,*) 'Adding more active blocks:'
          endif

          do while (nblocks_active < tasks)

             block_is_active_new(:) = block_is_active(:)

             do nb = 0, nblocks-1
                if (.not. block_is_active(nb)) then
                   if (block_is_active(west_block(nb))  .or. block_is_active(east_block(nb)) .or.  &
                       block_is_active(south_block(nb)) .or. block_is_active(north_block(nb))) then
                      if (verbose_active_blocks) write(iulog,*) 'Activate block', nb
                      block_is_active_new(nb) = .true.
                      nblocks_active = nblocks_active + 1
                      if (nblocks_active == tasks) exit
                   endif
                endif
             enddo   ! nb

             block_is_active(:) = block_is_active_new(:)

          enddo   ! nblocks_active < tasks

          if (nblocks_active < tasks) then  ! should not happen, but check just in case
             write(iulog,*) 'Error, still have nblocks_active < tasks:', nblocks_active, tasks
             call parallel_stop(__FILE__, __LINE__)
          endif

          if (verbose_active_blocks) then
             write(iulog,*) 'After activating more blocks: nblocks, nblocks_active:', nblocks, nblocks_active
             write(iulog,*) ' '
             write(iulog,*) 'Block layout:'
             do j = nstasks-1, 0, -1
                do i = 0, ewtasks-1
                   nb = ewtasks*j + i
                   write(iulog, '(i5)', advance='no') nb
                enddo
                write(iulog,*) ' '
             enddo
             write(iulog,*) ' '
             write(iulog,*) 'block_is_active:'
             do j = nstasks-1, 0, -1
                do i = 0, ewtasks-1
                   nb = ewtasks*j + i
                   write(iulog, '(l5)', advance='no') block_is_active(nb)
                enddo
                write(iulog,*) ' '
             enddo
          endif  ! verbose_active_blocks

       endif   ! nblocks_active < tasks

    endif   ! main_task

    ! Now that we have one task per active block, set up grid info on the local task.

    ! Broadcast active blocks to all tasks
    call broadcast(block_is_active)
    nblocks_active = count(block_is_active)

    ! Set up a correspondence between tasks and active blocks
    allocate(block_to_task(0:nblocks-1))
    allocate(task_to_block(0:tasks-1))

    block_to_task(:) = -1   ! do not set to 0, since 0 is a valid block and task number
    task_to_block(:) = -1

    nt = -1
    do nb = 0, nblocks-1
       if (block_is_active(nb)) then
          nt = nt + 1
          task_to_block(nt) = nb
          block_to_task(nb) = nt
       endif
    enddo

    if (main_task .and. verbose_active_blocks) then
       write(iulog,*) ' '
       write(iulog,*) 'Task layout:'
       do j = nstasks-1, 0, -1
          do i = 0, ewtasks-1
             nb = ewtasks*j + i
             write(iulog, '(i5)', advance='no') block_to_task(nb)
          enddo
          write(iulog,*) ' '
       enddo
    endif

    ! Assign grid info for the local task (nt = this_rank)
    ! This info can be copied directly from the block arrays computed above.
    nt = this_rank
    nb = task_to_block(nt)

    ewrank = ewrank_block(nb)
    global_col_offset = global_col_offset_block(nb)
    ewlb = ewlb_block(nb)
    ewub = ewub_block(nb)
    own_ewn = own_ewn_block(nb)
    local_ewn = local_ewn_block(nb)
    ewn = ewn_block(nb)

    nsrank = nsrank_block(nb)
    global_row_offset = global_row_offset_block(nb)
    nslb = nslb_block(nb)
    nsub = nsub_block(nb)
    own_nsn = own_nsn_block(nb)
    local_nsn = local_nsn_block(nb)
    nsn = nsn_block(nb)

    ! Assign W, E, S and N neighbor tasks, based on the block layout.
    ! Note: The halo logic for no_ice BCs is designed so that if (say) west = this_rank,
    !  then no halo values are sent to the west neighbor.

    if (block_is_active(west_block(nb))) then   ! set 'west' to the task that owns west_block(nb)
       west = block_to_task(west_block(nb))
    else  ! set 'west' to the local task
       west = nt
    endif

    if (block_is_active(east_block(nb))) then   ! set 'east' to the task that owns east_block(nb)
       east = block_to_task(east_block(nb))
    else  ! set 'east' to the local task
       east = nt
    endif

    if (block_is_active(south_block(nb))) then   ! set 'south' to the task that owns south_block(nb)
       south = block_to_task(south_block(nb))
    else  ! set 'south' to the local task
       south = nt
    endif

    if (block_is_active(north_block(nb))) then   ! set 'north' to the task that owns north_block(nb)
       north = block_to_task(north_block(nb))
    else  ! set 'north' to the local task
       north = nt
    endif

    if (main_task .and. verbose_active_blocks) then
       write(iulog,*) ' '
       write(iulog,*) 'this_rank, nb, west(nb), east(nb), south(nb), north(nb):',  &
                this_rank, nb, west_block(nb), east_block(nb), south_block(nb), north_block(nb)
       write(iulog,*) 'this_rank, nt, west(nt), east(nt), south(nt), north(nt):',  &
                this_rank, nt, west, east, south, north
    endif

    ! Identify blocks at the corner of the global domain.  Each such block has a locally owned
    ! corner cell where scalars should be set to zero in halo updates.

    southwest_corner = .false.
    southeast_corner = .false.
    northeast_corner = .false.
    northwest_corner = .false.

    if (this_rank > west .and. this_rank > south) then  ! west and south blocks are active
       corner_block = west_block(south_block(nb))
       if (.not.block_is_active(corner_block)) then     ! southwest block is not active
          southwest_corner = .true.
       endif
    endif

    if (this_rank < east .and. this_rank > south) then  ! east and south blocks are active
       corner_block = east_block(south_block(nb))
       if (.not.block_is_active(corner_block)) then     ! southeast block is not active
          southeast_corner = .true.
       endif
    endif

    if (this_rank < east .and. this_rank < north) then  ! east and north blocks are active
       corner_block = east_block(north_block(nb))
       if (.not.block_is_active(corner_block)) then     ! northeat block is not active
          northeast_corner = .true.
       endif
    endif

    if (this_rank > west .and. this_rank < north) then  ! west and north blocks are active
       corner_block = west_block(north_block(nb))
       if (.not.block_is_active(corner_block)) then     ! northwest block is not active
          northwest_corner = .true.
       endif
    endif

    if (verbose_active_blocks) then
!       if (southwest_corner) write(iulog,*) 'Southwest corner, task =', this_rank
!       if (southeast_corner) write(iulog,*) 'Southeast corner, task =', this_rank
!       if (northeast_corner) write(iulog,*) 'Northeast corner, task =', this_rank
!       if (northwest_corner) write(iulog,*) 'Northwest corner, task =', this_rank
    endif

    ! Set the limits of locally owned vertices on the staggered grid.
    ! For periodic BC, staggered_ilo = staggered_jlo = staggered_lhalo+1.
    ! For no_ice BC, we use the same values as for periodic.
    ! This means that vertices along the southern and western global boundaries
    !  are not considered to be locally owned by any task.  This is acceptable since
    !  there is no ice adjacent to the boundary, and thus (u,v) = 0 for boundary vertices.

    staggered_ilo = staggered_lhalo+1
    staggered_ihi = ewn - 1 - staggered_uhalo
    staggered_jlo = staggered_lhalo+1
    staggered_jhi = nsn - 1 - staggered_uhalo

    ! Check that we have not split up the problem too much.  We do not want halos overlapping in either dimension.
    ! local_* - lhalo - uhalo is the actual number of non-halo cells on a processor.
    if ((local_nsn - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
        write(iulog,*) "NS halos overlap on processor ", this_rank
        call parallel_stop(__FILE__, __LINE__)
    endif

    if ((local_ewn  - lhalo - uhalo) .lt. (lhalo + uhalo + 1)) then
        write(iulog,*) "EW halos overlap on processor ", this_rank
        call parallel_stop(__FILE__, __LINE__)
    endif

    if (verbose_active_blocks) then
    call parallel_barrier
!       write(iulog,*) 'task, west, east, south, north:', this_rank, west, east, south, north
!       write(iulog,*) 'task, SW, SE, NE, NW:', this_rank, &
!            southwest_corner, southeast_corner, northwest_corner, northeast_corner
    endif

    if (present(reprosum_in)) then
       reprosum = reprosum_in
    else
       reprosum = .false.
    endif

    ! If computing reproducible sums, then set some options
    !TODO - Are these saved from one call to the next?
    ! Note: For standalone CISM, reprosum = F by default; can set = T in the config file
    !       For CESM coupled runs, reprosum = T by default

    if (reprosum) then

#ifdef CCSM_COUPLED
!!       call shr_reprosum_setops()
#else
!!       call cism_reprosum_setops()
#endif

    endif  ! reprosum

! Uncomment to print grid geometry
!    write(iulog,*) " "
!    write(iulog,*) "Process ", this_rank, " Total = ", tasks, " ewtasks = ", ewtasks, " nstasks = ", nstasks
!    write(iulog,*) "Process ", this_rank, " ewrank = ", ewrank, " nsrank = ", nsrank
!    write(iulog,*) "Process ", this_rank, " l_ewn = ", local_ewn, " o_ewn = ", own_ewn
!    write(iulog,*) "Process ", this_rank, " l_nsn = ", local_nsn, " o_nsn = ", own_nsn
!    write(iulog,*) "Process ", this_rank, " ewlb = ", ewlb, " ewub = ", ewub
!    write(iulog,*) "Process ", this_rank, " nslb = ", nslb, " nsub = ", nsub
!    write(iulog,*) "Process ", this_rank, " east = ", east, " west = ", west
!    write(iulog,*) "Process ", this_rank, " north = ", north, " south = ", south
!    write(iulog,*) "Process ", this_rank, " ew_vars = ", own_ewn, " ns_vars = ", own_nsn
!    write(iulog,*) "Process ", this_rank, " global_col_offset = ", global_col_offset, &
!                                      " global_row_offset = ", global_row_offset

    call distributed_print_grid(own_ewn, own_nsn)

    end associate

  end subroutine distributed_grid_active_blocks

!=======================================================================

  subroutine distributed_print_grid(l_ewn,l_nsn)

    ! Gathers and prints the overall grid layout by processor counts.
    use mpi_mod
    implicit none

    integer :: l_ewn, l_nsn

    integer :: i,j,curr_count
    integer,dimension(2) :: mybounds
    integer,dimension(:,:),allocatable :: bounds

    ! begin
    mybounds(1) = l_ewn
    mybounds(2) = l_nsn

    if (main_task) then
       allocate(bounds(2,tasks))
    else
       allocate(bounds(1,1))
    end if
    call fc_gather_int(mybounds,2,mpi_integer,bounds,2,mpi_integer,main_rank,comm)
    if (main_task) then
       do i = 1,tasks
          if (bounds(1,i) .ne. -1) then
             ! total up number of processors with matching distribution
             curr_count = 1
             do j = i+1,tasks
                if ((bounds(1,i) .eq. bounds(1,j)) .and. (bounds(2,i) .eq. bounds(2,j))) then
                   ! if matching current distribution, increment counter
                   curr_count = curr_count + 1
                   bounds(1,j) = -1  ! mark so not counted later
                   bounds(2,j) = -1
                endif
             enddo
             write(iulog,*) "Layout(EW,NS) = ", bounds(1,i), bounds(2,i), " total procs = ", curr_count
          endif
       end do
    end if
    ! automatic deallocation

  end subroutine distributed_print_grid

!=======================================================================

  !WHL - Added the next two subroutines to create communicators for rows and columns of tasks.
  !      These are useful for a tridiagonal solve in 1D along a global row or column.
  !      Cf. parallel_initialise and parallel_set_info

  subroutine parallel_create_comm_row(comm, parallel)

    use mpi_mod
    implicit none
    integer, intent(in) :: comm          ! global communicator
    type(parallel_type) :: parallel

    integer, parameter :: my_main_rank_row = 0
    integer :: ierror

    associate(  &
         nsrank        => parallel%nsrank,        &
         comm_row      => parallel%comm_row,      &
         tasks_row     => parallel%tasks_row,     &
         this_rank_row => parallel%this_rank_row, &
         main_rank_row => parallel%main_rank_row, &
         main_task_row => parallel%main_task_row  &
         )

    ! Create a communicator (comm_row) that groups rows according to their value of nsrank,
    !  and assign a main_task to each row.

    call mpi_comm_split(comm, nsrank, this_rank, comm_row, ierror)
    call mpi_comm_size(comm_row, tasks_row, ierror)
    call mpi_comm_rank(comm_row, this_rank_row, ierror)
    main_rank_row = my_main_rank_row
    main_task_row = (this_rank_row==main_rank_row)

    if (main_task_row) write(iulog,*) 'Create comm_row: this_rank, tasks_row, this_rank_row, main_task_row =', &
         this_rank, tasks_row, this_rank_row, main_task_row

    end associate

  end subroutine parallel_create_comm_row

!=======================================================================

  subroutine parallel_create_comm_col(comm, parallel)

    use mpi_mod
    implicit none
    integer, intent(in) :: comm          ! global communicator
    type(parallel_type) :: parallel

    integer, parameter :: my_main_rank_col = 0
    integer :: ierror

    associate(  &
         ewrank        => parallel%ewrank,        &
         comm_col      => parallel%comm_col,      &
         tasks_col     => parallel%tasks_col,     &
         this_rank_col => parallel%this_rank_col, &
         main_rank_col => parallel%main_rank_col, &
         main_task_col => parallel%main_task_col  &
         )

    ! Create a communicator (comm_col) that groups columns according to their value of ewrank,
    !  and assign a main_task to each column.

    call mpi_comm_split(comm, ewrank, this_rank, comm_col, ierror)
    call mpi_comm_size(comm_col, tasks_col, ierror)
    call mpi_comm_rank(comm_col, this_rank_col, ierror)
    main_rank_col = my_main_rank_col
    main_task_col = (this_rank_col==main_rank_col)

    if (main_task_col) write(iulog,*) 'Create comm_col: this_rank, tasks_col, this_rank_col, main_task_col =', &
         this_rank, tasks_col, this_rank_col, main_task_col

    end associate

  end subroutine parallel_create_comm_col

!=======================================================================

  subroutine parallel_localindex(iglobal, jglobal, ilocal, jlocal, rlocal, parallel)

    ! Calculates the local i,j indices and rank from the global i,j indices
    integer,intent(in) :: iglobal, jglobal 
    integer,intent(out)  :: ilocal, jlocal, rlocal
    type(parallel_type) :: parallel

    integer :: flag, flag_out

    associate(  &
         own_ewn           => parallel%own_ewn,            &
         own_nsn           => parallel%own_nsn,            &
         global_row_offset => parallel%global_row_offset,  &
         global_col_offset => parallel%global_col_offset   &
         )

    flag = 0   ! This flag will be flipped on exactly one processor if the global point is valid
    ilocal = iglobal + lhalo - global_col_offset
    jlocal = jglobal + lhalo - global_row_offset

    ! Check whether these are valid values of ilocal and jlocal
    ! If so, then flip the flag and broadcast these values
    if ( (ilocal > lhalo .and. ilocal <= lhalo + own_ewn)  &
                         .and.                             &
         (jlocal > lhalo .and. jlocal <= lhalo + own_nsn) ) then
       flag = 1
    endif

    call parallel_reduce_maxloc(flag, flag_out, rlocal)

    if (flag_out==1) then
       call broadcast(ilocal, rlocal)
       call broadcast(jlocal, rlocal)
    else ! global indices are invalid
       if (main_task) then
          write(iulog,*) 'Invalid global indices: iglobal, jglobal =', iglobal, jglobal
          call parallel_stop(__FILE__,__LINE__)
       endif
    endif

    end associate

  end subroutine parallel_localindex

!=======================================================================

  subroutine parallel_test_comm_row_col(parallel)

    ! Test the communicators for rows and columns of processors.
    ! Optionally, these communicators can be used when solving a 1D system of equations
    !  along a single row or column of the domain.

    use mpi_mod

    type(parallel_type), intent(in) :: parallel

    integer :: i, j, ierror
    real(dp), dimension(:,:), allocatable :: test_array
    real(dp), dimension(:,:), allocatable :: global_test_array

    ! Test the row and column communicators

    ! row gather/scatter

    ! mpi barriers not strictly needed, but can make print statements cleaner
    call mpi_barrier(parallel%comm_row,ierror)
    call mpi_barrier(parallel%comm_col,ierror)

    allocate(test_array(2,parallel%own_nsn))

    do j = 1, parallel%own_nsn
       do i = 1, 2
          test_array(i,j) = (this_rank + 2) * real(i*j, dp)
       enddo
    enddo

    if (this_rank == 0 .or. this_rank == 1) then
       write(iulog,*) ' '
       do i = 1, 2
          write(iulog,*) 'test_array, i, this_rank =', i, this_rank
          do j = 1, parallel%own_nsn
             write(iulog,'(f6.0)',advance='no') test_array(i,j)
          enddo
          write(iulog,*) ' '
          write(iulog,*) ' '
       enddo
    endif   ! this_rank

    call gather_var_row(test_array, global_test_array, parallel)

!!    if (parallel%main_task_row) then
    if (parallel%main_task_row .and. this_rank == 0) then
       do i = 1, 2*parallel%tasks_row
          write(iulog,*) 'Gathered row test_array, this_rank, i =', this_rank, i
          do j = 1, size(global_test_array,2)
             write(iulog,'(f6.0)',advance='no') global_test_array(i,j)
          enddo
          write(iulog,*) ' '
       enddo
       write(iulog,*) ' '
    endif

    call scatter_var_row(test_array, global_test_array, parallel)

    if (this_rank == 0) then
       write(iulog,*) ' '
       do i = 1, 2
          write(iulog,*) 'Scattered test_array, i, this_rank =', i, this_rank
          do j = 1, parallel%own_nsn
             write(iulog,'(f6.0)',advance='no') test_array(i,j)
          enddo
          write(iulog,*) ' '
          write(iulog,*) ' '
       enddo
    endif   ! this_rank

    deallocate(test_array)

    ! column gather/scatter

    call mpi_barrier(parallel%comm_row,ierror)
    call mpi_barrier(parallel%comm_col,ierror)

    allocate(test_array(2,parallel%own_ewn))
    do j = 1, 2
       do i = 1, parallel%own_ewn
          test_array(j,i) = (this_rank + 2) * real(j*i, dp)
       enddo
    enddo

    if (this_rank == 0 .or. this_rank == 2) then
       write(iulog,*) ' '
       do j = 1, 2
          write(iulog,*) 'test_array: j, this_rank =', j, this_rank
          do i = 1, parallel%own_ewn
             write(iulog,'(f6.0)',advance='no') test_array(j,i)
          enddo
          write(iulog,*) ' '
          write(iulog,*) ' '
       enddo
    endif   ! this_rank

    call gather_var_col(test_array, global_test_array, parallel)

!!    if (parallel%main_task_col) then
    if (parallel%main_task_col .and. this_rank == 0) then
       do j = 1, 2*parallel%tasks_col
          write(iulog,*) 'Gathered column test_array, this_rank, j =', this_rank, j
          do i = 1, size(global_test_array,2)
             write(iulog,'(f6.0)',advance='no') global_test_array(j,i)
          enddo
          write(iulog,*) ' '
       enddo
       write(iulog,*) ' '
    endif

    call scatter_var_col(test_array, global_test_array, parallel)

    if (this_rank == 0) then
       write(iulog,*) ' '
       do j = 1, 2
          write(iulog,*) 'Scattered test_array, j, this_rank =', j, this_rank
          do i = 1, parallel%own_ewn
             write(iulog,'(f6.0)',advance='no') test_array(j,i)
          enddo
          write(iulog,*) ' '
          write(iulog,*) ' '
       enddo
    endif   ! this_rank

    deallocate(test_array)

  end subroutine parallel_test_comm_row_col

!=======================================================================

end module cism_parallel_decomp

!=======================================================================
