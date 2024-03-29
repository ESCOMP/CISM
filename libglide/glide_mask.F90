!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_mask.F90 - part of the Community Ice Sheet Model (CISM)  
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

#include "glide_mask.inc"

module glide_mask

    ! masking ice thicknesses

    use glimmer_global, only : dp
    use cism_parallel, only: lhalo, uhalo, parallel_type, parallel_halo, parallel_reduce_sum

    implicit none

contains

!TODO - Remove iarea and ivol calculations?  They are now computed in glide_write_diag..

!TODO - Write a new subroutine (in addition to glide_set_mask) to compute mask for staggered grid?
!       This subroutine is now called from glissade_velo_driver with stagthck and stagtopg
!       as input arguments.

  subroutine glide_set_mask(numerics, thck, topg, ewn, nsn, eus, mask, iarea, ivol, &
                            exec_serial, parallel)

    use glide_types
    use glimmer_physcon, only : rhoi, rhoo
    implicit none

    type(glide_numerics), intent(in) :: numerics !Numerical parameters structure
    real(dp), dimension(:,:), intent(in) :: thck ! Ice thickness
    real(dp), dimension(:,:), intent(in) :: topg ! Bedrock topography (not lower surface!)
    integer, intent(in) :: ewn, nsn              ! Grid size
    real(dp), intent(in) :: eus                  ! Sea level
    integer, dimension(:,:), intent(inout) :: mask   ! Output mask
    real(dp), intent(inout), optional :: ivol, iarea ! Area and volume of ice

    logical, optional :: exec_serial             !JEFF If executing in serial in MPI program.
    type(parallel_type), optional :: parallel    ! info for parallel communication

    ! local variables
    integer ew,ns
    logical :: exec_serial_flag

    !Note - This array may not be needed, at least in parallel.

    ! Create an array to "fake" the boundaries of the mask so that boundary
    ! finding can work even on the boundaries of the real mask.

    integer, dimension(0:ewn+1,0:nsn+1) :: maskWithBounds;

    !TODO - What is the exec_serial option?  Is it still needed?

    !JEFF Handle exec_serial optional parameter
    if ( present(exec_serial) ) then
       exec_serial_flag = exec_serial
    else
       ! Default to off
       exec_serial_flag = .FALSE.
    endif

    mask = 0

    if (present(iarea)) iarea = 0.d0
    if (present(ivol)) ivol = 0.d0

!Note - This mask is confusing.  Wondering if we should replace it by a series of logical masks.

! Would need the following:
! glide_mask_has_ice = 1
! glide_mask_thin_ice = 3
! glide_mask_ocean = 4 (below sea level, with or without ice)
! glide_mask_land = 8 (complement of glide_mask_ocean)
! glide_mask_grounding_line = 16 (could define in terms of margin and has ice?)
! glide_mask_margin = 32 (has_ice + at least one neighbor with no ice)
! glide_mask_dirichlet_bc = 64
! glide_mask_comp_domain_bnd = 128 (no longer needed with new global BC?)
! glide_no_ice (complement of glide_has_ice)
! glide_is_thin
! glide_is_ocean (ocean + no_ice; change to glide_ocean_icefree or remove?)
! glide_is_land (land + no_ice; change to glide_land_icefree or remove?)
! glide_is_ground (land + has_ice)
! glide_is_float (ocean + has_ice)
! glide_is_grounding_line (just inside or just outside? Used only in glide_ground)
! glide_is_margin
! glide_is_land_margin (margin + land + has_ice)
! glide_is_calving (margin + ocean + has_ice; change the name to is_marine_margin?)
! glide_is_marine_ice_edge (margin + (float or GL); may not be needed)
! glide_is_dirichlet_boundary 
! glide_is_comp_domain_bnd (may not be needed with new global BC?)
! 
!       If we keep the present structure, could change glide_is_land to glide_icefree_land,
!                                                      glide_is_ocean to glide_icefree_ocean
!       Could get by with fewer masks in the code by removing some combinations
!       Could remove *BITS

    !Identify points with any ice
    where (thck > 0.d0)
        mask = ior(mask, GLIDE_MASK_HAS_ICE)  ! GLIDE_MASK_HAS_ICE = 1; see glide_mask.inc
    endwhere

    !Identify points where the ice is below the ice dynamics limit
    where (thck > 0.d0 .and. thck < numerics%thklim)
        mask = ior(mask, GLIDE_MASK_THIN_ICE)  ! GLIDE_MASK_THIN_ICE = 3
    endwhere

    !Identify points where the ice is floating or where there is open ocean
    where (topg - eus < (-rhoi/rhoo) * thck)
        mask = ior(mask, GLIDE_MASK_OCEAN)   ! GLIDE_MASK_OCEAN = 8
    elsewhere
        mask = ior(mask, GLIDE_MASK_LAND)    ! GLIDE_MASK_LAND = 4
    endwhere

    if (present(iarea) .and. present(ivol)) then
        call get_area_vol(thck, numerics%dew, numerics%dns, numerics%thklim, iarea, ivol, exec_serial_flag)
    end if

    !TODO - Replace the following with a halo call for 'mask', with appropriate global BC?
    
    maskWithBounds = 0
    maskWithBounds(1:ewn, 1:nsn) = MASK
    maskWithBounds(0,1:nsn) = mask(1,:)
    maskWithBounds(1:ewn,0) = mask(:,1)
    maskWithBounds(ewn+1,1:nsn) = mask(ewn,:)
    maskWithBounds(1:ewn,nsn+1) = mask(:,nsn)
    maskWithBounds(0,0) = mask(1,1)
    maskWithBounds(ewn+1,nsn+1) = mask(ewn,nsn)
    maskWithBounds(0,nsn+1) = mask(1,nsn)
    maskWithBounds(ewn+1,0) = mask(ewn,1)

    ! finding boundaries

    !Note: If halo cells are present, maskWithBounds array may not be needed; can replace with mask array.
    !      Not sure what happens here when we're computing a mask on the velocity grid.

    do ns = 1,nsn
       do ew = 1,ewn
          !Find the grounding line
          if (GLIDE_IS_GROUND(MASK(ew,ns))) then    ! land + has_ice
             if (GLIDE_IS_FLOAT(maskWithBounds(ew-1,ns)) .or. &
                  GLIDE_IS_FLOAT(maskWithBounds(ew+1,ns)) .or. &
                  GLIDE_IS_FLOAT(maskWithBounds(ew,ns-1)) .or. & 
                  GLIDE_IS_FLOAT(maskWithBounds(ew,ns+1))) then
                MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_GROUNDING_LINE)
             end if
          end if

          ! Ice margin
          ! *tb* A point is now masked even if it touches the ocean on one corner.
          if ( GLIDE_HAS_ICE(mask(ew, ns)) .and. &
              (GLIDE_NO_ICE(maskWithBounds(ew-1,ns))   .or. GLIDE_NO_ICE(maskWithBounds(ew+1,ns))   .or. &
               GLIDE_NO_ICE(maskWithBounds(ew,ns-1))   .or. GLIDE_NO_ICE(maskWithBounds(ew,ns+1))   .or. &
               GLIDE_NO_ICE(maskWithBounds(ew-1,ns-1)) .or. GLIDE_NO_ICE(maskWithBounds(ew-1,ns+1)) .or. &
               GLIDE_NO_ICE(maskWithBounds(ew+1,ns-1)) .or. GLIDE_NO_ICE(maskWithBounds(ew+1,ns+1)))) then
             MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_MARGIN)
          end if

       end do
    end do

    !JEFF Don't call halo update if running in serial mode
    !WHL - I think the halo update will now work in serial mode.
    if (.NOT. exec_serial_flag .and. present(parallel)) then
       call parallel_halo(mask, parallel)
    endif

  end subroutine glide_set_mask

  subroutine augment_kinbc_mask(mask, kinbcmask)

    !  Augments the Glide mask with the location of kinematic (dirichlet) boundary
    !   conditions.  These locations cannot be determined by the model a priori, and
    !   must be specified through a field in a NetCDF file.
    integer, dimension(:,:), target :: mask
    integer, dimension(:,:) :: kinbcmask

    integer, dimension(:,:), pointer :: maskp

    !Because the kinematic boundary conditions are specified on the staggered grid,
    !there may be a size mismatch here depending on whether we are computing a mask
    !for the staggered grid.
    if (size(mask, 1) /= size(kinbcmask, 1)) then
        maskp => mask(1:size(mask,1) - 1, 1:size(mask,2) - 1)
    else
        maskp => mask
    end if

    where (kinbcmask /= 0)
        maskp = ior(maskp, GLIDE_MASK_DIRICHLET_BC)
    endwhere
  end subroutine augment_kinbc_mask

  subroutine get_area_vol(thck, dew, dns, thklim, iarea, ivol, exec_serial)

    use glimmer_paramets, only : len0, thk0

    implicit none
    real(dp), dimension(:,:) :: thck
    real(dp) :: dew, dns, thklim
    real(dp) :: iarea, ivol, sum(2)
    logical :: exec_serial

    integer :: i,j

    do i = 1+lhalo, size(thck,1)-uhalo
        do j = 1+lhalo, size(thck,2)-uhalo
            if (thck(i,j) > thklim ) then
                iarea = iarea + 1
                ivol = ivol + thck(i,j)
            end if
        end do
    end do

    iarea = iarea  * dew * dns
    ivol = ivol * dew * dns
    
    if (.NOT. exec_serial) then
       sum(1) = iarea
       sum(2) = ivol
       sum = parallel_reduce_sum(sum)
       iarea = sum(1)
       ivol  = sum(2)
    endif

    ! convert from model units to SI units
    iarea = iarea*len0*len0
    ivol = ivol*len0*len0*thk0

  end subroutine get_area_vol
 
  subroutine calc_iareaf_iareag(dew, dns, mask, iareaf, iareag, exec_serial)
    
    use glimmer_paramets, only : len0 
    implicit none
    real(dp), intent(in) :: dew, dns
    real(dp), intent(out) :: iareaf, iareag
    integer, dimension(:,:), intent(in) :: mask 
    logical, optional :: exec_serial  ! If executing in serial in MPI program.

    integer :: i,j
    logical :: exec_serial_flag
    real(dp) :: sum(2)
 
    !TODO - exec_serial option may not be needed
    if ( present(exec_serial) ) then
      exec_serial_flag = exec_serial
    else
      ! Default to off
      exec_serial_flag = .FALSE.
    endif

    iareaf = 0.d0
    iareag = 0.d0 

    !loop over locally owned scalars
    do j = 1+lhalo, size(mask,2)-uhalo
      do i = 1+lhalo, size(mask,1)-uhalo
        if (GLIDE_IS_FLOAT(mask(i,j))) then
          iareaf = iareaf + dew * dns
        else if(GLIDE_IS_GROUND_OR_GNDLINE(mask(i,j))) then
          iareag = iareag + dew * dns
        end if
      end do
    end do

    if (.NOT. exec_serial_flag) then
       sum(1) = iareaf
       sum(2) = iareag
       sum = parallel_reduce_sum(sum)
       iareaf = sum(1)
       iareag = sum(2)
    endif

    ! convert from model units to SI units (m^2)
    iareag = iareag*len0*len0
    iareaf = iareaf*len0*len0

  end subroutine calc_iareaf_iareag

    subroutine glide_marine_margin_normal(thck, mask, marine_bc_normal, &
                                          exec_serial, parallel)

      !TODO - Remove subroutine glide_marine_margin_normal?  Old PBJ routine.
      !       Also can remove calc_normal_45deg

        use glimmer_physcon, only:pi
        implicit none
        !> This subroutine derives from the given mask the normal to an ice shelf
        !> each point on the marine margin.
        real(dp), dimension(:,:), intent(in) :: thck
        integer, dimension(:,:), intent(in) :: mask
        real(dp), dimension(:,:), intent(out) :: marine_bc_normal

        logical, optional :: exec_serial  !JEFF If executing in serial in MPI program
        type(parallel_type), optional :: parallel    ! info for parallel communication

        integer :: i, j, dx, dy, k
        logical :: exec_serial_flag

        real(dp), dimension(size(thck,1), size(thck,2)) :: direction_x, direction_y
        
        real(dp), dimension(-1:1, -1:1) :: angle_lookup

        !JEFF Handle exec_serial optional parameter
        if ( present(exec_serial) ) then
           exec_serial_flag = exec_serial
        else
           ! Default to off
           exec_serial_flag = .FALSE.
        endif

        !direction_y =    -1       0       1        !direction_x = 
        angle_lookup(-1, :) = (/ 3*pi/4,   pi/2,   pi/4 /)  !-1
        angle_lookup( 0, :) = (/   pi,     0D0,  2*pi   /)  ! 0
        angle_lookup( 1, :) = (/ 5*pi/4, 3*pi/2, 7*pi/4 /)  ! 1
        call upwind_from_mask(mask, direction_x, direction_y, exec_serial_flag)

        !Set up a thickness variable with "ghost cells" so that we don't go out
        !of bounds with the vectorized operation below
        !thckWithBounds(1:size(thck,1), 1:size(thck,2)) = thck
        !thckWithBounds(:,0) = thckWithBounds(:,1)
        !thckWithBounds(0,:) = thckWithBounds(1,:)
        !thckWithBounds(size(thck,1)+1,:) = thckWithBounds(size(thck,1),:)
        !thckWithBounds(:,size(thck,2)+1) = thckWithBounds(:,size(thck,2))
        do i = 1, size(mask, 1)
            do j = 1, size(mask, 2)
                if (GLIDE_IS_CALVING(mask(i,j))) then
                    dx = int(direction_x(i,j))
                    dy = int(direction_y(i,j))
                    if (dx == 0 .and. dy == 0) then
                        write(*,*)"A shelf front point has been identified at:"
                        write(*,*)"x = ",i
                        write(*,*)"y = ",j
                        write(*,*)"But neither x nor y derivatives have been marked as upwinded."
                        write(*,*)"This should never happen, if this error appears it is a bug"
                        write(*,*)"and should be reported."
                        write(*,*)"The mask around this point follows:"
                        write(*,*)"--------------------------"

                        !Write a header row with a * in the column corresponding to the center
                        do k = -4, 4
                            if (k==0) then
                                write(*,"(A)",advance="no")"           *"
                            else if (i+k > 0 .and. i+k <= size(mask,1)) then
                                write(*,"(A)",advance="no")"            "
                            end if
                        end do 
                        write(*,*)

                        do k=4, -4, -1
                            if (j+k > 0 .and. j+k <= size(mask, 2)) then
                                if (k == 0) then
                                    write(*,*) "*", mask(max(1,i-4):min(size(mask,1),i+4),j+k)
                                else
                                    write(*,*) " ", mask(max(1,i-4):min(size(mask,1),i+4),j+k)
                                end  if
                            end if
                        end do
                        write(*,*)"--------------------------"
                        write(*,*)"Have a nice day!"
                        !stop
                    end if
                    marine_bc_normal(i,j) = angle_lookup(dx, dy) 
                    !marine_bc_normal(i,j) = calc_normal_45deg(thckWithBounds(i-1:i+1,j-1:j+1))
                else
                   ! NOTE(wjs, 2016-01-07) this was set to NaN, but I am setting it to 0
                   ! because I'm removing nan_mod. Since this routine is unused, this
                   ! seemed to be an okay thing to do.
                   marine_bc_normal(i,j) = 0._dp
                end if
            end do
        end do
        if (.NOT. exec_serial_flag .and. present(parallel)) then
           call parallel_halo(marine_bc_normal, parallel)
        endif
      end subroutine glide_marine_margin_normal

    function calc_normal_45deg(thck3x3)
        use glimmer_physcon, only: pi
        
        !> Computes the angle of the normal vector, in radians, for the given
        !> 3x3 segment of ice geometry.
        !> The normal is given in increments of 45 degrees (no nicer
        !> interpolation is currently done)
        !> This is based on the Payne and Price GLAM code, if/when this is
        !> integrated into CISM it should probably be refactored to use this.
        real(dp), dimension(3,3) :: thck3x3

        real(dp) :: calc_normal_45deg
         
        real(dp), dimension(3,3) :: mask, maskcorners
        real(dp), dimension(3,3) :: thckmask
        real(dp), dimension(3) :: testvect
        real(dp) :: phi, deg2rad
        integer :: loc_latbc

        deg2rad = pi / 180.0d0
        loc_latbc = 0
        phi = 0.d0
        mask(:,1) = (/ 0.0d0, 180.0d0, 0.0d0 /)
        mask(:,2) = (/ 270.0d0, 0.0d0, 90.0d0 /)
        mask(:,3) = (/ 0.0d0, 360.0d0, 0.0d0 /)
        maskcorners(:,1) = (/ 225.0d0, 0.0d0, 135.0d0 /)
        maskcorners(:,2) = (/ 0.0d0, 0.0d0, 0.0d0 /)
        maskcorners(:,3) = (/ 315.0d0, 0.0d0, 45.0d0 /)

        ! specify new value of 'loc' vector such that fwd/bwd diffs. are set up correctly in sparse matrix
        ! when function 'fillsprsebndy' is called. Also, specify appropriate values for the vectors 'normal'
        ! and 'fwdorbwd', which specify the orientation of the boundary normal and the direction of forward or
        ! backward differencing to be done in the lateral boundary condition functions 'normhorizmainbc_lat'
        ! and 'crosshorizmainbc_lat'

        ! following is algorithm for calculating boundary normal at 45 deg. increments, based on arbitray
        ! boundary shape

        where( thck3x3  /=  0.0d0 )
            thckmask = 0.0_dp
        elsewhere( thck3x3 == 0.0d0 )
            thckmask = 1.0d0
        endwhere

        testvect = sum( thckmask * mask, 1 )

        !if( up == 3 )then ! temporary code for debugging
        !  do i = 3,1,-1
        !  print *, 'thck = ', thck(:,i)
        !  end do
        !  print *, ' '
        !
        !  do i = 3,1,-1
        !      print *, 'thckmask = ', thckmask(:,i)
        !  end do
        !  print *, ' '
        !
        !  print *, 'testvect =  ', testvect
        !  print *, ' '
        !end if

        ! calculate the angle of the normal in cart. (x,y) system w/ 0 deg. at 12 O'clock, 90 deg. at 3 O'clock, etc.
        if( sum( sum( thckmask, 1 ) ) == 1.0d0 )then
            phi = sum( sum( thckmask * maskcorners, 1 ) )
        else
            if( any( testvect == 360.0d0 ) )then
                if( sum( testvect ) == 450.0d0 )then
                    phi = 45.0d0
                elseif( sum( testvect ) == 630.0d0 )then
                    phi = 315.0d0
                else
                    phi = 0.0d0
                end if
            elseif( all( testvect  /=  360 ) )then
                phi = sum( testvect ) / sum( testvect/testvect, testvect  /=  0.0d0 )
            end if
        end if

        calc_normal_45deg = deg2rad * phi
        
        !Tim's Note: This appears to actually compute 0 at 6 O'clock according
        !to Glimmer's coordinate system.  90 deg. is still 3 O'clock.
        !I'm going to correct for this here rather than dig through the code
        !above
        calc_normal_45deg = pi - calc_normal_45deg 
        if (calc_normal_45deg < 0) calc_normal_45deg = calc_normal_45deg + 2*pi

    end function

!TODO - Remove subroutine upwind_from_mask?  Not currently used.

    !Fills a field of differencing directions suitable to give a field
    !derivative routine.  Uses centered differencing everywhere except for the
    !marine ice margin, where upwinding and downwinding is used to avoid
    !differencing across the boundary.

    subroutine upwind_from_mask(mask, direction_x, direction_y, &
                                exec_serial, parallel)

        integer, dimension(:,:), intent(in) :: mask
        double precision, dimension(:,:), intent(out) :: direction_x, direction_y

        logical, optional :: exec_serial  !JEFF If executing in serial in MPI program.
        type(parallel_type), optional :: parallel    ! info for parallel communication

        integer :: i,j
        logical :: exec_serial_flag

        !JEFF Handle exec_serial optional parameter
        if ( present(exec_serial) ) then
           exec_serial_flag = exec_serial
        else
           ! Default to off
           exec_serial_flag = .FALSE.
        endif

        direction_x = 0
        direction_y = 0

        !Detect locations of the marine margin
        do i = 1, size(mask,1)
            do j = 1, size(mask,2)
                if (GLIDE_IS_CALVING(mask(i,j))) then
                    !Detect whether we need to upwind or downwind in the Y
                    !direction
                    if (i > 1) then
                        if (.not. GLIDE_HAS_ICE(mask(i-1,j))) then
                            direction_x(i,j) = 1
                        end if
                    end if

                    if (i < size(mask, 1)) then
                        if (.not. GLIDE_HAS_ICE(mask(i+1,j))) then
                            direction_x(i,j) = -1
                        end if
                    end if

                    !Detect whether we need to upwind or downwind in the X
                    !direction
                    if (j > 1) then
                        if (.not. GLIDE_HAS_ICE(mask(i,j-1))) then
                            direction_y(i,j) = 1
                        end if
                    end if
                    
                    if (j < size(mask, 2)) then
                        if (.not. GLIDE_HAS_ICE(mask(i,j+1))) then
                            direction_y(i,j) = -1
                        end if
                    end if

                    !If we are at a point that is "interior" to two other boundary points, 
                    !such as the lower right of:
                    !o b i
                    !b b i
                    !(o = ocean, b = boundary, i = interior), then we will not detect the need
                    !to upwind or downwind.  However, we still should for consistency with other
                    !mask points (in some cases, not doing so can lead to a singular calculation
                    !at the marine ice front)
                    !
                    !We can think of this operation as avoiding calving points where there is 
                    !a non-calving point to upwind into.
                    !
                    !NOTE: We need a better way to detect interior points.  Right now I am just using
                    !points that are floating, and that works, but this doesn't work for two reasons:
                    !1. Boundary points are also floating
                    !2. Could fail for a very thin ice shelf
                    if (int(direction_x(i,j)) == 0 .and. int(direction_y(i,j)) == 0 .and. &
                        i > 1 .and. j > 1 .and. i < size(mask, 1) .and. j < size(mask, 2)) then
                        if (.not. GLIDE_HAS_ICE(mask(i-1, j-1))) then
                            direction_x(i,j) = 1
                            direction_y(i,j) = 1
                        else if (.not. GLIDE_HAS_ICE(mask(i-1, j+1))) then
                            direction_x(i,j) = 1
                            direction_y(i,j) = -1
                        else if (.not. GLIDE_HAS_ICE(mask(i+1, j-1))) then
                            direction_x(i,j) = -1
                            direction_y(i,j) = 1
                        else if (.not. GLIDE_HAS_ICE(mask(i+1, j+1))) then
                            direction_x(i,j) = -1
                            direction_y(i,j) = -1
                        end if
                    end if
                end if
            end do
        end do

        if (.NOT. exec_serial_flag .and. present(parallel)) then
            call parallel_halo(direction_x, parallel)
            call parallel_halo(direction_y, parallel)
        endif

    end subroutine upwind_from_mask

end module glide_mask
