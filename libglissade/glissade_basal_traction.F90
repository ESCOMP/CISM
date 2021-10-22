!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_basal_traction.F90 - part of the Community Ice Sheet Model (CISM)  
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

#include "glide_mask.inc"
#include "config.inc"

  module glissade_basal_traction

  !-----------------------------------------------------------------------------
  ! Compute or prescribe the basal traction coefficient 'beta' as required by
  ! the higher-order velocity solver.
  ! 
  ! Note that beta is assumed to be a positive constant.  In earlier versions of
  ! the code it was called 'betasquared'.
  !
  ! The units are Pa/(m/yr) if we assume a linear sliding law of the form
  !    taub_x = -beta * u, taub_y = -beta * v
  !
  ! However, the units are Pa if beta is treated as a till yield stress.
  !
  ! See glide_types.F90 for the meaning of the various options.
  !
  !-----------------------------------------------------------------------------

  use glimmer_paramets, only : dp
  use glimmer_physcon,  only : scyr
  use glimmer_paramets, only : vel0, tau0
  use glimmer_log
  use glide_types
  use cism_parallel, only : this_rank, main_task, parallel_type, &
       parallel_halo, staggered_parallel_halo, parallel_globalindex, distributed_scatter_var

  implicit none

  private
  public :: calcbeta, calc_effective_pressure, glissade_init_effective_pressure

!***********************************************************************

contains

!***********************************************************************

  subroutine calcbeta (whichbabc,                    &
                       parallel,                     &
                       dew,           dns,           &
                       ewn,           nsn,           &
                       thisvel,       othervel,      &
                       basal_physics,                &
                       flwa_basal,    thck,          &
                       topg,          eus,           &
                       ice_mask,                     &
                       land_mask,                    &
                       f_ground,                     &
                       beta_external,                &
                       beta,                         &
                       which_ho_beta_limit,          &
                       which_ho_powerlaw_c,          &
                       which_ho_coulomb_c,           &
                       itest, jtest,  rtest)

  ! subroutine to calculate map of beta sliding parameter, based on 
  ! user input ("whichbabc" flag, from config file as "which_ho_babc").
   
  ! NOTE: Previously, the input arguments were assumed to be dimensionless
  ! and were rescaled in this routine.  Now the input arguments are
  ! assumed to have the units given below.
     
  use glimmer_paramets, only: len0
  use glimmer_physcon, only: gn, pi
  use glissade_grid_operators, only: glissade_stagger

  implicit none

  ! Input/output arguments

  integer, intent(in) :: whichbabc     ! option for basal friction parameterization

  type(parallel_type), intent(in) :: &
       parallel                        ! info for parallel communication

  integer, intent(in) :: ewn, nsn      ! horizonal grid dimensions

  real(dp), intent(in)                    :: dew, dns           ! m
  real(dp), intent(in), dimension(:,:)    :: thisvel, othervel  ! basal velocity components (m/yr)
  type(glide_basal_physics), intent(in)   :: basal_physics      ! basal physics object
  real(dp), intent(in), dimension(:,:)    :: flwa_basal         ! flwa for the basal ice layer (Pa^{-3} yr^{-1})
  real(dp), intent(in), dimension(:,:)    :: thck               ! ice thickness (m)
  real(dp), intent(in), dimension(:,:)    :: topg               ! bed topography (m)
  real(dp), intent(in)                    :: eus                ! eustatic sea level (m) relative to z = 0

  integer, intent(in), dimension(:,:) :: &
       ice_mask,        & ! = 1 where ice is present (thck > thklim), else = 0
       land_mask          ! = 1 where topg > eus

  real(dp), intent(in), dimension(:,:)    :: f_ground           ! grounded ice fraction at vertices, 0 <= f_ground <= 1
  real(dp), intent(in), dimension(:,:)    :: beta_external      ! fixed beta read from external file (Pa yr/m)
  real(dp), intent(inout), dimension(:,:) :: beta               ! basal traction coefficient (Pa yr/m)

  integer, intent(in) :: which_ho_beta_limit                    ! option to limit beta for grounded ice
                                                                ! 0 = absolute based on beta_grounded_min; 1 = weighted by f_ground
  integer, intent(in) :: which_ho_powerlaw_c                    ! basal friction option for Cp
  integer, intent(in) :: which_ho_coulomb_c                     ! basal frection option for Cc
  integer, intent(in), optional :: itest, jtest, rtest          ! coordinates of diagnostic point

  ! Local variables

  ! Note: Adding fields for parallel ISHOM-C test case
  real(dp), dimension(:,:), allocatable :: beta_global          ! beta on the global grid
  real(dp), dimension(:,:), allocatable :: beta_extend          ! beta extended to the ice grid (dimensions ewn, nsn)

  real(dp) :: smallnum = 1.0d-2  ! m/yr

  real(dp) :: Ldomain   ! size of full domain
  real(dp) :: omega     ! frequency of beta field
  real(dp) :: dx, dy
  integer :: ew, ns, i, j

  real(dp), dimension(size(beta,1), size(beta,2)) :: speed      ! ice speed, sqrt(uvel^2 + vvel^2), m/yr

  ! variables for power law
  real(dp) :: powerlaw_p, powerlaw_q

  ! variables for Coulomb friction law
  real(dp) :: coulomb_c   ! Coulomb law friction coefficient (unitless)
  real(dp) :: powerlaw_c_const  ! power law friction coefficient (Pa m^{-1/3} yr^{1/3})
  real(dp) :: lambda_max        ! wavelength of bedrock bumps at subgrid scale (m)
  real(dp) :: m_max             ! maximum bed obstacle slope (unitless)
  real(dp) :: m                 ! exponent m in power law

  integer, dimension(size(thck,1), size(thck,2)) :: &
       ice_or_land_mask,   & ! = 1 where ice_mask = 1 or land_mask = 1, else = 0       
       imask                 ! = 1 where thck > 0, else = 1

  real(dp), dimension(size(beta,1), size(beta,2)) ::  &
       big_lambda,         & ! bedrock characteristics
       flwa_basal_stag       ! basal flwa interpolated to the staggered grid (Pa^{-n} yr^{-1})

  ! variables for Tsai et al. parameterization
  real(dp) :: taub_powerlaw  ! basal shear stress given by a power law as in Tsai et al. (2015)
  real(dp) :: taub_coulomb   ! basal shear stress given by Coulomb friction as in Tsai et al. (2015)

  ! variables for pseudo-plastic law
  real(dp) :: q              ! exponent for pseudo-plastic law (unitless)
                             ! q = 1 for linear, q = 0 for plastic, 0 < q < 1 for power law
  real(dp) :: u0             ! threshold velocity for pseudo-plastic law (m/yr)
  real(dp) :: phi            ! phi for pseudoplastic law (degress, 0 < phi < 90)
  real(dp) :: tanphi         ! tan(phi) for pseudo-plastic law (unitless)
  real(dp) :: bed            ! bed elevation, topg - eus (m)
  real(dp) :: phimin, phimax ! min and max values of phi for pseudo-plastic law (degrees)
  real(dp) :: bedmin, bedmax ! bed elevations (m) below which phi = phimin and above which phi = phimax
  real(dp) :: tau_c          ! yield stress for pseudo-plastic law (unitless)
  real(dp) :: numerator, denominator

  character(len=300) :: message

  integer :: iglobal, jglobal

  real(dp) :: effecpress_capped   ! capped effective pressure for Coulomb laws (ZI specifically)

  logical, parameter :: verbose_beta = .false.

  ! Compute the ice speed: used in power laws where beta = beta(u).
  ! Enforce a minimum speed to prevent beta from become very large when velocity is small.
  speed(:,:) = dsqrt(thisvel(:,:)**2 + othervel(:,:)**2 + smallnum**2)

  ! If beta_powerlaw_umax is set to a nonzero value, then limit the speed to this value.
  ! Note: The actual ice speed can be greater than umax.  This is just a way of shutting off the feedback
  !        between beta and ice speed (beta down as speed up) when the ice speed is large.
  !       It helps make the model more stable.
  if (basal_physics%beta_powerlaw_umax > 0.0d0) then
     speed(:,:) = min(speed(:,:), basal_physics%beta_powerlaw_umax)
  endif

  ! Compute coulomb_c; used in basal friction laws with yield stress proportional to coulomb_c

  if (which_ho_coulomb_c == HO_COULOMB_C_CONSTANT) then
     ! set coulomb_c = constant value
     basal_physics%coulomb_c(:,:) = basal_physics%coulomb_c_const
  elseif (which_ho_coulomb_c == HO_COULOMB_C_ELEVATION) then

     ! set coulomb_c based on bed elevation
     call set_coulomb_c_elevation(ewn,        nsn,   &
                                  topg,       eus,   &
                                  basal_physics,     &
                                  basal_physics%coulomb_c)

  else  ! HO_COULOMB_C_INVERSION, HO_COULOMB_C_EXTERNAL
     ! do nothing; use coulomb_c as computed elsewhere
  endif

  ! Compute powerlaw_c; used in basal friction laws with beta proportional to u^(1/m)

  if (which_ho_powerlaw_c == HO_POWERLAW_C_CONSTANT) then
     ! set powerlaw_c = constant value
     basal_physics%powerlaw_c(:,:) = basal_physics%powerlaw_c_const
  else  ! HO_POWERLAW_C_INVERSION, HO_POWERLAW_C_EXTERNAL
     ! do nothing; use powerlaw_c as computed elsewhere
  endif

  ! Compute beta based on whichbabc

  select case(whichbabc)

    case(HO_BABC_BETA_CONSTANT)   ! spatially uniform beta value; useful for debugging and test cases

          beta(:,:) = basal_physics%ho_beta_const  ! Pa yr/m

    case(HO_BABC_BETA_BPMP)  ! large value for frozen bed, lower value for bed at pressure melting point

          where(basal_physics%bpmp_mask == 1)      ! bed is at pressure melting point
             beta(:,:) = basal_physics%ho_beta_small    ! Pa yr/m
          elsewhere 
             beta(:,:) = basal_physics%ho_beta_large    ! Pa yr/m
          endwhere

    case(HO_BABC_PSEUDO_PLASTIC)  ! pseudo-plastic sliding law using the new coulomb_c options

       ! Pseudo-plastic sliding law from PISM:
       !
       ! (tau_bx,tau_by) = -tau_c * (u,v) / (u_0^q * |u|^(1-q))
       ! where the yield stress tau_c = N * tan(phi), or equivalently tau_c = N * coulomb_c
       ! N = effective pressure, computed in subroutine calc_effective_pressure
       ! q, u0 and phi are user-configurable parameters:
       !    q = exponent (q = 1 for linear sliding, q = 0 for a plastic bed, 0 < q < 1 for power-law behavior), default = 1/3
       !    u0 = threshold velocity (the velocity at which tau_b = tau_c), default = 100 m/yr
       !    0 < coulomb_c < 1
       ! As in PISM, coulomb_c is allowed to vary with bed elevation.
       ! See Aschwanden et al. (2013), The Cryosphere, 7, 1083-1093, Supplement; see also the PISM Users Guide.

       q  = basal_physics%pseudo_plastic_q
       u0 = basal_physics%pseudo_plastic_u0

       ! compute beta based on N, coulomb_c and u
       do ns = 1, nsn-1
          do ew = 1, ewn-1
             tau_c = basal_physics%effecpress_stag(ew,ns) * basal_physics%coulomb_c(ew,ns)
             beta(ew,ns) = tau_c / (u0**q * speed(ew,ns)**(1.0d0 - q))

             !WHL - debug
             if (verbose_beta .and. present(rtest) .and. present(itest) .and. present(jtest)) then
                if (this_rank == rtest .and. ew == itest .and. ns == jtest) then
                   write(6,*) 'i, j, bed, coulomb_c, tau_c, speed, beta:', &
                        ew, ns, bed, phi, basal_physics%coulomb_c(ew,ns), tau_c, speed(ew,ns), beta(ew,ns)
                endif
             endif
          enddo   ! ew
       enddo   ! ns

    case(HO_BABC_PSEUDO_PLASTIC_OLD)  ! older method, retained for backward compatibility
                                      ! TODO: Remove the old method when no longer the CESM default

       q  = basal_physics%pseudo_plastic_q
       u0 = basal_physics%pseudo_plastic_u0

       phimin = basal_physics%pseudo_plastic_phimin
       phimax = basal_physics%pseudo_plastic_phimax
       bedmin = basal_physics%pseudo_plastic_bedmin
       bedmax = basal_physics%pseudo_plastic_bedmax

       ! Note: There is a minor bug in the loop below.
       !       beta(ew,ns) is computed based on topg(ew,ns); should use stagtopg(ew,ns) instead.
       !       Leaving the bug as is for back compatibility, given that this method will be deprecated.

       do ns = 1, nsn-1
          do ew = 1, ewn-1

             ! compute tan(phi) based on bed elevation
             bed = topg(ew,ns) - eus
             if (bed <= bedmin) then
                phi = phimin
             elseif (bed >= bedmax) then
                phi = phimax
             else   ! bed elevation is between bedmin and bedmax
                phi = phimin + ((bed - bedmin)/(bedmax - bedmin)) * (phimax - phimin)
             endif
             tanphi = tan(phi * pi/180.d0)

             ! compute beta based on tan(phi), N and u
             tau_c = tanphi * basal_physics%effecpress_stag(ew,ns) 
             beta(ew,ns) = tau_c / (u0**q * speed(ew,ns)**(1.0d0 - q))

             !WHL - debug
             if (verbose_beta .and. present(rtest) .and. present(itest) .and. present(jtest)) then
                if (this_rank == rtest .and. ew == itest .and. ns == jtest) then
                   write(6,*) 'i, j, bed, tanphi, tau_c, speed, beta:', &
                        ew, ns, bed, tanphi, tau_c, speed(ew,ns), beta(ew,ns)
                endif
             endif

          enddo   ! ew
       enddo   ! ns

    case(HO_BABC_YIELD_PICARD)  ! take input value for till yield stress and force beta to be implemented such
                                ! that plastic-till sliding behavior is enforced (see additional notes in documentation).

      !!! NOTE: Eventually, this option could provide the till yield stress as calculate from the basal processes submodel.
      !!!       Currently, to enable sliding over plastic till, simply specify the value of "beta" as 
      !!!       if it were the till yield stress (in units of Pascals).
      
      beta(:,:) = basal_physics%mintauf(:,:)*tau0 &                                      ! plastic yield stress (converted to Pa)
                         / dsqrt( thisvel(:,:)**2 + othervel(:,:)**2 + (smallnum)**2 )   ! velocity components (m/yr)

      !!! since beta is updated here, communicate that info to halos
      call staggered_parallel_halo(beta, parallel)

    case(HO_BABC_BETA_LARGE)      ! frozen (u=v=0) ice-bed interface

       !Note: This option is redundant in that it could be implemented using HO_BETA_CONST
       !      But keeping it for historical reasons since many config files use it

       beta(:,:) = basal_physics%ho_beta_large      ! Pa yr/m  (= 1.0d10 by default)

    case(HO_BABC_ZOET_IVERSON)

       ! Use the sliding law proposed by Zoet & Iverson (2020):
       !     tau_b = N * C_c * [u_b/(u_b + u_t)]^(1/m), Eq. 3 in ZI(2020)
       ! where N   = effective pressure
       !       C_c = a constant in the range [0,1]
       !       u_t = threshold speed controlling the transition between powerlaw and Coulomb behavior
       !       m   = powerlaw exponent
       !Note: We have added the option to cap N at a value of N_max.
       !      By default, N_max is large enough that there will be no limiting,
       !       but N_max can be set to a smaller value in the config file..

       m = basal_physics%powerlaw_m

       do ns = 1, nsn-1
          do ew = 1, ewn-1
             effecpress_capped = min(basal_physics%effecpress_stag(ew,ns), &
                                     basal_physics%zoet_iverson_nmax)
             tau_c = basal_physics%coulomb_c(ew,ns) * effecpress_capped
             beta(ew,ns) = tau_c * speed(ew,ns)**(1.0d0/m - 1.0d0)  &
                  / (speed(ew,ns) + basal_physics%zoet_iverson_ut)**(1.0d0/m)

             !WHL - debug
             if (verbose_beta .and. present(rtest) .and. present(itest) .and. present(jtest) .and. &
                  this_rank == rtest .and. ew == itest .and. ns == jtest) then
                write(6,*) 'Cc, N, speed, beta =', basal_physics%coulomb_c(ew,ns), &
                     basal_physics%effecpress_stag(ew,ns), speed(ew,ns), beta(ew,ns)
             endif

          enddo
       enddo

    case(HO_BABC_ISHOMC)          ! prescribe according to ISMIP-HOM test C

       !TODO: Carry out this operation at initialization, before calling calcbeta?
       !Note: Ideally, beta would be read in from an external netCDF file.
       !      However, this is not possible given that the global velocity grid is smaller
       !       than the ice grid and hence not able to fit the full beta field.
       !      The following code sets beta on the full grid as prescribed by Pattyn et al. (2008).

       ! Allocate a global array on the main task only.
       ! On other tasks, allocate a size 0 array, since distributed_scatter_var wants to deallocate on all tasks.
       if (main_task) then
          allocate(beta_global(parallel%global_ewn, parallel%global_nsn))
       else
          allocate(beta_global(0,0))
       endif

       ! Prescribe beta as in Pattyn et al., The Cryosphere, 2008
       ! Note: These beta values live at vertices, not cell centers.
       !       They need a global array of size (ewn,nsn) to hold values on the global boundary.
       if (main_task) then

          Ldomain = parallel%global_ewn * dew   ! size of full domain (must be square)
          omega = 2.d0*pi / Ldomain

          beta_global(:,:) = 0.0d0
          do ns = 1, parallel%global_nsn
             do ew = 1, parallel%global_ewn
                dx = dew * ew
                dy = dns * ns
                beta_global(ew,ns) = 1000.d0 + 1000.d0 * sin(omega*dx) * sin(omega*dy)
             enddo
          enddo

       endif

       ! Scatter the global beta values back to local arrays
       ! Note: beta_extend has dimensions (ewn,nsn), so it can receive scattered data from beta_global.
       allocate(beta_extend(ewn, nsn))
       beta_extend(:,:) = 0.d0
       call distributed_scatter_var(beta_extend, beta_global, parallel)

       ! distributed_scatter_var does not update the halo, so do an update here
       call parallel_halo(beta_extend, parallel)

       ! Copy beta_extend to beta on the local processor.
       ! This is done since beta lives on the velocity grid and has dimensions (ewn-1,nsn-1).
       beta(:,:) = 0.d0
       do ns = 1, nsn-1
          do ew = 1, ewn-1
             beta(ew,ns) = beta_extend(ew, ns)
          enddo
       enddo

      ! beta_extend is no longer needed (beta_global is deallocated in distributed_scatter_var)
      deallocate(beta_extend)

    case(HO_BABC_BETA_EXTERNAL)   ! use beta value from external file

       ! set beta to the prescribed external value
       ! Note: This assumes that beta_external has units of Pa yr/m on input.
       beta(:,:) = beta_external(:,:)

       ! beta is initialized to a negative value; we can use that fact to check whether
       ! it has been read correctly from the file
       if (maxval(beta) < 0.d0) then
          call write_log('ERROR: Trying to use HO_BABC_EXTERNAL_BETA, but all beta values are < 0,')
          call write_log('which implies that beta could not be read from the input file.')
          call write_log('Make sure that beta is in the cism input file,')
          call write_log('or change which_ho_babc to a different option.')
          call write_log('Invalid value for beta. See log file for details.', GM_FATAL)
       end if

    case(HO_BABC_POWERLAW)   ! a simple power law
       !   Assume taub = C * ub^(1/m)
       ! implying beta = C * ub^(1/m - 1) 
       ! m should be a positive exponent

       do ns = 1, nsn-1
          do ew = 1, ewn-1
             beta(ew,ns) = basal_physics%powerlaw_c(ew,ns) &
                         * speed(ew,ns)**(1.0d0/basal_physics%powerlaw_m - 1.0d0)

             !WHL - debug
             if (verbose_beta .and. present(rtest) .and. present(itest) .and. present(jtest)) then
                if (this_rank == rtest .and. ew == itest .and. ns == jtest) then
                   write(6,*) 'r, i, j, Cp, speed, beta:', &
                        rtest, itest, jtest, basal_physics%powerlaw_c(ew,ns), speed(ew,ns), beta(ew,ns)
                endif
             endif
          enddo
       enddo

    case(HO_BABC_POWERLAW_EFFECPRESS)   ! a power law that uses effective pressure
       !TODO - Remove POWERLAW_EFFECPRESS option? Rarely if ever used.
       ! See Cuffey & Paterson, Physics of Glaciers, 4th Ed. (2010), p. 240, eq. 7.17
       ! This is based on Weertman's classic sliding relation (1957) augmented by the bed-separation index described by Bindschadler (1983)
       !   ub = k taub^p N^-q
       ! rearranging for taub gives:
       !   taub = k^(-1/p) ub^(1/p) N^(q/p)

       ! p and q should be _positive_ exponents. If p/=1, this is nonlinear in velocity.
       ! Cuffey & Paterson recommend p=3 and q=1, and k dependent on thermal & mechanical properties of ice and inversely on bed roughness.   
       !TODO - Change powerlaw_p to powerlaw_m, and make powerlaw_q a config parameter

       powerlaw_p = 3.0d0
       powerlaw_q = 1.0d0

       beta(:,:) = basal_physics%friction_powerlaw_k**(-1.0d0/powerlaw_p) &
            * basal_physics%effecpress_stag(:,:)**(powerlaw_q/powerlaw_p) &
            * speed(:,:)**(1.0d0/powerlaw_p - 1.0d0)

    case(HO_BABC_COULOMB_FRICTION)

      ! TODO: Remove this option; effectively the same as the Schoof option below
      !       Might need to modify MISMIP test config files that use this option

      ! Basal stress representation using Schoof sliding law with Coulomb friction
      ! See Schoof 2005 PRS, eqn. 6.2  (see also Pimentel, Flowers & Schoof 2010 JGR)

       ! Set up parameters needed for the friction law
       m_max = basal_physics%coulomb_bump_max_slope       ! maximum bed obstacle slope(unitless)
       lambda_max = basal_physics%coulomb_bump_wavelength ! wavelength of bedrock bumps (m)

       ! Need flwa of the basal layer on the staggered grid
       !TODO - Pass in ice_mask instead of computing imask here?
       !       (Small difference: ice_mask = 1 where thck > thklim rather than thck > 0)
       where (thck > 0.d0)
          imask = 1
       elsewhere
          imask = 0
       end where
       call glissade_stagger(ewn,         nsn,               &
                             flwa_basal,  flwa_basal_stag,   &
                             imask,       stagger_margin_in = 1)
       ! TODO Not sure if a halo update is needed on flwa_basal_stag!  I don't think so if nhalo>=2.

       ! Compute biglambda = wavelength of bedrock bumps [m] * flwa [Pa^-n yr^-1] / max bed obstacle slope [dimensionless]
       big_lambda(:,:) = (lambda_max / m_max) * flwa_basal_stag(:,:)

       ! Note: For MISMIP3D, coulomb_c is multiplied by a spatial factor (c_space_factor) which is
       !       read in during initialization. This factor is typically between 0 and 1.
       !       If this factor is not present in the input file, it is set to 1 everywhere.

       ! Compute beta
       ! Note: Where this equation has powerlaw_m, we used to have Glen's flow exponent n,
       !       following the notation of Leguy et al. (2014).
       !       Changed to powerlaw_m to be consistent with the Schoof and Tsai laws.
       m = basal_physics%powerlaw_m
       beta(:,:) = basal_physics%coulomb_c(:,:) * basal_physics%effecpress_stag(:,:) &
            * speed(:,:)**(1.0d0/m - 1.0d0) * &
            (speed(:,:) + basal_physics%effecpress_stag(:,:)**m * big_lambda)**(-1.0d0/m)

       ! If c_space_factor /= 1.0 everywhere, then multiply beta by c_space_factor
       ! TODO: Replace c_space_factor with a spatially varying coulomb_c field.
       if (maxval(abs(basal_physics%c_space_factor_stag(:,:) - 1.0d0)) > tiny(0.0d0)) then
          beta(:,:) = beta(:,:) * basal_physics%c_space_factor_stag(:,:)
       endif

       ! Limit for numerical stability
       !TODO - Is limiting needed?
       where (beta > 1.0d8)
          beta = 1.0d8
       end where

    case(HO_BABC_COULOMB_POWERLAW_SCHOOF)

       ! Use the basal friction formulation of Schoof (2005), modified following Asay-Davis et al. (2016).
       ! This formulation uses a constant value of basal flwa, which allows several Coulomb parameters
       !  (lambda_max, m_max and flwa_basal) to be combined into a single parameter powerlaw_c, 
       !  as in the Tsai power law below.
       !
       ! The equation for tau_b = beta * u_b is
       ! 
       !                    powerlaw_c * coulomb_c * N
       ! tau_b = ---------------------------------------------- u_b^{1/m}
       !         [powerlaw_c^m * u_b + (coulomb_c * N)^m]^{1/m}
       !
       ! where m = powerlaw_m
       !
       ! This is the second modified basal traction law in MISMIP+. See Eq. 11 of Asay-Davis et al. (2016).
       ! Note: powerlaw_c corresponds to beta^2 in their notation, and coulomb_c corresponds to alpha^2.
       !
       ! Depending on the value of which_ho_powerlaw_c and which_ho_coulomb_c, there are different ways
       !  to apply this sliding law:
       ! (0) Set powerlaw_c and coulomb_c to a constant everywhere.
       ! (1) Obtain spatially varying powerlaw_c or coulomb_c fields by inversion.
       ! (2) Use spatially varying powerlaw_c or coulomb_c fields prescribed from a previous inversion.
       !
       ! Note: This law and the Tsai law are often run with spatially varying powerlaw_c,
       !  but have not yet been tested with spatially varying coulomb_c.

       m = basal_physics%powerlaw_m

       do ns = 1, nsn-1
          do ew = 1, ewn-1

             numerator = basal_physics%powerlaw_c(ew,ns) * basal_physics%coulomb_c(ew,ns)  &
                       * basal_physics%effecpress_stag(ew,ns)
             denominator = (basal_physics%powerlaw_c(ew,ns)**m * speed(ew,ns) +  &
                  (basal_physics%coulomb_c(ew,ns) * basal_physics%effecpress_stag(ew,ns))**m )**(1.d0/m)
             beta(ew,ns) = (numerator/denominator) * speed(ew,ns)**(1.d0/m - 1.d0)

             !WHL - debug
             if (verbose_beta .and. present(rtest) .and. present(itest) .and. present(jtest)) then
                if (this_rank == rtest .and. ew == itest .and. ns == jtest) then
                   print*, ' '
                   write(6,*) 'r, i, j, Cp, denom_u, denom_N, speed, beta, taub:', &
                        rtest, ew, ns, basal_physics%powerlaw_c(ew,ns), &
                        (basal_physics%powerlaw_c(ew,ns)**m * speed(ew,ns))**(1.d0/m), &
                        (basal_physics%coulomb_c(ew,ns) * basal_physics%effecpress_stag(ew,ns)), &
                        speed(ew,ns), beta(ew,ns), beta(ew,ns)*speed(ew,ns)
                endif
             endif

          enddo
       enddo

       ! If c_space_factor /= 1.0 everywhere, then multiply beta by c_space_factor
       if (maxval(abs(basal_physics%c_space_factor_stag(:,:) - 1.0d0)) > tiny(0.0d0)) then
          beta(:,:) = beta(:,:) * basal_physics%c_space_factor_stag(:,:)
       endif

       ! Limit for numerical stability
       !TODO - Is limiting needed?
       where (beta > 1.0d8)
          beta = 1.0d8
       end where

       !WHL - debug - Write values along a flowline
!      write(6,*) ' '
!      write(6,*) 'Apply Coulomb friction: i, j, speed, N_stag, beta, taub:'
!      ns = jtest
!      do ew = itest, itest+15
!         write(6,*) ew, ns, speed(ew,ns), basal_physics%effecpress_stag(ew,ns), beta(ew,ns), beta(ew,ns)*speed(ew,ns)
!      enddo

    case(HO_BABC_COULOMB_POWERLAW_TSAI)

      ! Basal stress representation based on Tsai et al. (2015)
      ! The basal stress is the minimum of two values:
      ! (1) power law:          tau_b = powerlaw_c * |u_b|^(1/powerlaw_m)
      ! (2) Coulomb friction:   tau_b = coulomb_c * N
      !                             N = effective pressure = rhoi*g*(H - H_f)
      !                           H_f = flotation thickness = (rhow/rhoi)*(eus-topg)
      ! This value of N is obtained by setting p_ocean_penetration = 1.0 in the config file.
      ! The other parameters (powerlaw_c, powerlaw_m and coulomb_c) can also be set in the config file.

       do ns = 1, nsn-1
          do ew = 1, ewn-1
             
             taub_powerlaw = basal_physics%powerlaw_c(ew,ns) * speed(ew,ns)**(1.d0/basal_physics%powerlaw_m)
             taub_coulomb  = basal_physics%coulomb_c(ew,ns) * basal_physics%effecpress_stag(ew,ns)

             if (taub_coulomb <= taub_powerlaw) then   ! apply Coulomb stress, which is smaller
                beta(ew,ns) = taub_coulomb / speed(ew,ns)
             else  ! apply power-law stress
                beta(ew,ns) = taub_powerlaw / speed(ew,ns)
             endif

          enddo   ! ew
       enddo   ! ns

       ! If c_space_factor /= 1.0 everywhere, then multiply beta by c_space_factor
       if (maxval(abs(basal_physics%c_space_factor_stag(:,:) - 1.0d0)) > tiny(0.0d0)) then
          beta(:,:) = beta(:,:) * basal_physics%c_space_factor_stag(:,:)
       endif

    case(HO_BABC_SIMPLE)    ! simple pattern; also useful for debugging and test cases
                            ! (here, a strip of weak bed surrounded by stronger bed to simulate an ice stream)

      beta(:,:) = 1.d4        ! Pa yr/m

      !TODO - Change this loop to work in parallel (set beta on the global grid and scatter to local)
      do ns = 5, nsn-5
         do ew = 1, ewn-1
            beta(ew,ns) = 100.d0      ! Pa yr/m
         end do
      end do

    case default
       ! do nothing

   end select


   ! Multiply beta by f_ground to reduce friction at partly grounded vertices.
   ! Note: With a GLP, f_ground will have values between 0 and 1 at vertices adjacent to the GL.
   !       Without a GLP, f_ground = 0 or 1 everywhere based on a flotation criterion.
   !       By convention, f_ground = 1 for land, and f_ground = 0 for ice-free ocean.
   !
   ! For beta close to 0 beneath grounded ice, it is possible to generate unrealistically fast flow.
   ! To prevent this, set beta to a minimum value beneath grounded ice.
   ! The default value of beta_grounded_min = 0.0, but can be set to a nonzero value in the config file.

   ! The limiting can be done either before or after multiplying by f_ground.
   ! If done after, then beta >= beta_grounded_min at all vertices, including lightly grounded vertices.
   !  This is the more stable option for settings where there can be large driving stresses near the GL.
   ! If done before, then beta -> 0 as f_ground -> 0, even if beta_grounded_min > 0.

   do ns = 1, nsn-1
      do ew = 1, ewn-1

         if (which_ho_beta_limit == HO_BETA_LIMIT_ABSOLUTE) then  ! absolute limit based on beta_grounded_min

            ! multiplication by f_ground before limiting

            beta(ew,ns) = beta(ew,ns) * f_ground(ew,ns)

            if (f_ground(ew,ns) > 0.d0 .and. beta(ew,ns) < basal_physics%beta_grounded_min) then
               beta(ew,ns) = basal_physics%beta_grounded_min
            endif

         elseif (which_ho_beta_limit == HO_BETA_LIMIT_FLOATING_FRAC) then  ! weighted by f_ground after limiting

            ! multiplication by f_ground after limiting

            if (f_ground(ew,ns) > 0.d0 .and. beta(ew,ns) < basal_physics%beta_grounded_min) then
               beta(ew,ns) = basal_physics%beta_grounded_min
            endif

            beta(ew,ns) = beta(ew,ns) * f_ground(ew,ns)

         endif   ! which_ho_beta_limit

      enddo
   enddo

   ! Bug check: Make sure beta >= 0
   ! This check will find negative values as well as NaNs
   do ns = 1, nsn-1
      do ew = 1, ewn-1 
         if (beta(ew,ns) >= 0.d0) then
            ! do nothing
         else
            call parallel_globalindex(ew, ns, iglobal, jglobal, parallel)
            write(message,*) 'Invalid beta value in calcbeta: this_rank, i, j, iglobal, jglobal, beta, f_ground:', &
                 this_rank, ew, ns, iglobal, jglobal, beta(ew,ns), f_ground(ew,ns)
            call write_log(trim(message), GM_FATAL)
         endif
      end do
   end do
   
   ! halo update
   !TODO - Move this halo update to a higher level?
   call staggered_parallel_halo(beta, parallel)

   !WHL - debug
   if (verbose_beta .and. present(rtest) .and. present(itest) .and. present(jtest)) then
      if (this_rank == rtest) then
         ew = itest; ns = jtest
         write(6,*) 'End of calcbeta, r, i, j, speed, f_ground, beta:', &
              rtest, ew, ns, speed(ew,ns), f_ground(ew,ns), beta(ew,ns)
      endif
   endif

  end subroutine calcbeta

!***********************************************************************

  subroutine glissade_init_effective_pressure(which_effecpress, basal_physics)

    ! Initialize calculations related to effective pressure.
    ! Currently, the only thing to do is initialize an array for
    !  option which_effecpress = HO_EFFECPRESS_BWATFLX.
    ! Note: f_effecpress should not be reset if restarting.
    !       Currently, this subroutine is called only when *not* restarting

    ! Input/output arguments

    integer, intent(in) :: &
         which_effecpress    ! input option for effective pressure

    type(glide_basal_physics), intent(inout) :: &
         basal_physics       ! basal physics object

    if (which_effecpress == HO_EFFECPRESS_BWATFLX) then
       basal_physics%f_effecpress(:,:) = 1.0d0
    endif

  end subroutine glissade_init_effective_pressure

!***********************************************************************

  subroutine calc_effective_pressure (which_effecpress,             &
                                      parallel,                     &
                                      ewn,           nsn,           &
                                      basal_physics, basal_hydro,   &
                                      ice_mask,      floating_mask, &
                                      thck,          topg,          &
                                      eus,                          &
                                      delta_bpmp,                   &
                                      bwat,          bwatflx,       &
                                      dt,                           &
                                      itest, jtest,  rtest)

    ! Calculate the effective pressure N at the bed.
    ! By default, N is equal to the overburden pressure, rhoi*g*H.
    ! Optionally, N can be reduced by the presence of water at the bed
    !  (btemp near bpmp, or nonzero bwat or bwatflx).
    ! N can also be reduced where there is a hydrological connection to the ocean,
    !  through weighting by (1 - Hf/H)^p (where Hf is the flotation thickness).

    use glimmer_physcon, only: rhoi, grav, rhoo
    use glissade_grid_operators, only: glissade_stagger

    implicit none

    ! Input/output arguments

    integer, intent(in) :: &
         which_effecpress    ! input option for effective pressure

    type(parallel_type), intent(in) :: &
         parallel            ! info for parallel communication

    integer, intent(in) :: &
         ewn, nsn            ! grid dimensions

    type(glide_basal_physics), intent(inout) :: &
         basal_physics       ! basal physics object
                             ! includes effecpress, effecpress_stag and various parameters

    type(glide_basal_hydro), intent(inout) :: &
         basal_hydro         ! basal hydro object
                             ! includes bwat and various parameters

    integer, dimension(:,:), intent(in) :: &
         ice_mask,         & ! = 1 where ice is present (thk > thklim), else = 0
         floating_mask       ! = 1 where ice is present and floating, else = 0
 
    !NOTE: If used, the following 2D fields (delta_bpmp, bwat, bwatflx, thck and topg) need to be correct in halos.

    real(dp), dimension(:,:), intent(in) ::  &
         thck,             & ! ice thickness (m)
         topg                ! bed topography (m)

    real(dp), intent(in) ::  &
         eus                 ! eustatic sea level (m) relative to z = 0

    real(dp), dimension(:,:), intent(in), optional ::  &
         delta_bpmp,       & ! Tpmp - T at the bed (K), used for HO_EFFECPRESS_BPMP option
         bwat,             & ! basal water thickness (m), used for HO_EFFECPRESS_BWAT option
         bwatflx             ! basal water flux at the bed (m/yr), used for HO_EFFECPRESS_BWATFLX option

    real(dp), intent(in), optional :: dt                     ! time step (yr)

    integer, intent(in), optional :: itest, jtest, rtest     ! coordinates of diagnostic point

    ! Local variables

    real(dp) :: &
         bpmp_factor,     &  ! factor between 0 and 1, used in linear ramp based on bpmp
         relative_bwat,   &  ! ratio bwat/bwat_threshold, limited to range [0,1]
         df_dt               ! rate of change of f_effecpress

    real(dp), dimension(ewn,nsn) ::  &
         overburden,            & ! overburden pressure, rhoi*g*H
         effecpress_ocean_p,    & ! pressure reduced by ocean connection
         f_pattyn_2d              ! rhoo*(eus-topg)/(rhoi*thck)
                                  ! = 1 at grounding line, < 1 for grounded ice, > 1 for floating ice

    real(dp) :: ocean_p           ! exponent in effective pressure parameterization, 0 <= ocean_p <= 1

    real(dp) :: f_pattyn          ! rhoo*(eus-topg)/(rhoi*thck)
    real(dp) :: f_pattyn_capped   ! f_pattyn capped to lie in range [0,1]

    real(dp) :: frac
    integer :: i, j

    logical, parameter :: verbose_effecpress = .false.

    ! Initialize the effective pressure N to the overburden pressure, rhoi*g*H

    overburden(:,:) = rhoi*grav*thck(:,:)
    basal_physics%effecpress(:,:) = overburden(:,:)

    select case(which_effecpress)

    case(HO_EFFECPRESS_OVERBURDEN)

       ! do nothing; already initialized to overburden

       ! Note: Here we assume (unrealistically) that N = rhoi*g*H even for floating ice.
       !       However, the basal friction coefficient (beta) will equal zero for floating ice
       !        since it is weighted by the grounded ice fraction.

    case(HO_EFFECPRESS_BPMP)

       if (present(delta_bpmp)) then

          ! Reduce N where the basal temperature is near the pressure melting point,
          !  as defined by delta_bpmp = bpmp - Tbed.
          ! N decreases from overburden for a frozen bed to a small value for a thawed bed.
          ! bpmp_factor = 0 where the bed is thawed (delta_bpmp <= 0)
          ! bpmp_factor = 1 where the bed is frozen (delta_bpmp >= effecpress_bpmp_threshold)
          ! 0 < bpmp_factor < 1 where 0 < delta_bpmp < bpmp_threshold

          do j = 1, nsn
             do i = 1, ewn

                bpmp_factor = max(0.0d0, min(1.0d0, delta_bpmp(i,j)/basal_physics%effecpress_bpmp_threshold))
                basal_physics%effecpress(i,j) = basal_physics%effecpress(i,j) * &
                     (basal_physics%effecpress_delta + bpmp_factor * (1.0d0 - basal_physics%effecpress_delta))

                !TODO - not sure this is needed, because of weighting by f_ground
                ! set to zero for floating ice
                if (floating_mask(i,j) == 1) basal_physics%effecpress(i,j) = 0.0d0

             enddo
          enddo

       endif   ! present(delta_bpmp)

    case(HO_EFFECPRESS_BWAT)

       if (present(bwat)) then

          ! Reduce N where basal water is present.
          ! N decreases from overburden for bwat = 0 to a small value for bwat = effecpress_bwat_threshold.

          do j = 1, nsn
             do i = 1, ewn
                if (bwat(i,j) > 0.0d0) then

                   relative_bwat = max(0.0d0, min(bwat(i,j)/basal_physics%effecpress_bwat_threshold, 1.0d0))

                   basal_physics%effecpress(i,j) = basal_physics%effecpress(i,j) * &
                        (basal_physics%effecpress_delta + &
                        (1.0d0 - relative_bwat) * (1.0d0 - basal_physics%effecpress_delta))

                end if
             enddo
          enddo

       endif   ! present(bwat)

       !TODO - Not needed?
       where (floating_mask == 1)
          ! set to zero for floating ice
          basal_physics%effecpress = 0.0d0
       end where

    case(HO_EFFECPRESS_BWATFLX)

       ! Note: The units of bwatflux are volume per unit area per unit time, i.e. m/yr.
       !       This is the rate at which bwat would increase if there were inflow but no outflow.

       if (present(bwatflx)) then

          ! Reduce N where the basal water flux is greater than zero.
          ! This is done by prognosing f_effecpress = effecpress/overburden:
          !          df/dt = [1 - f*(F/F0)] / tau
          !          where f = f_effecpress, F = bwatflx, F0 = effecpress_bwatflx_threshold,
          !              tau = effecpress_timescale
          ! The steady-state f < 1 when F > F0.
          ! As f decreases, the marginal effect of additional flux also decreases.

          do j = 1, nsn
             do i = 1, ewn
                if (bwatflx(i,j) > 0.0d0) then

                   df_dt = ( 1.0d0 - basal_physics%f_effecpress(i,j) * &
                        (bwatflx(i,j)/basal_physics%effecpress_bwatflx_threshold) ) / &
                        basal_physics%effecpress_timescale
                   basal_physics%f_effecpress(i,j) = basal_physics%f_effecpress(i,j) + df_dt * dt

                   ! Limit to be in the range [effecpress_delta, 1.0)
                   basal_physics%f_effecpress(i,j) = min(basal_physics%f_effecpress(i,j), 1.0d0)
                   basal_physics%f_effecpress(i,j) = max(basal_physics%f_effecpress(i,j), basal_physics%effecpress_delta)

                   basal_physics%effecpress(i,j) = basal_physics%f_effecpress(i,j) * overburden(i,j)

                end if
             enddo
          enddo

          if (verbose_effecpress .and. this_rank == rtest) then
             print*, ' '
             print*, 'After bwatflx, f_effecpress, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   if (thck(i,j) > 0.0d0) then
                      write(6,'(f10.5)',advance='no') basal_physics%f_effecpress(i,j)
                   else
                      write(6,'(f10.5)',advance='no') 0.0d0
                   endif
                enddo
                write(6,*) ' '
             enddo
          endif

       endif   ! present(bwatflx)

       !TODO - Modify for deluxe GLP?
       !TODO - Not sure this is needed, because beta is later weighted by f_ground.
!       where (floating_mask == 1)
!          ! set to zero for floating ice
!          basal_physics%effecpress = 0.0d0
!       end where

    case(HO_EFFECPRESS_BWAT_BVP)

       if (present(bwat)) then

          ! Reduce N where basal water is present, following Bueler % van Pelt (2015).
          ! N decreases from overburden P_0 for bwat = 0 to a small value for bwat = effecpress_bwat_threshold.
          ! This scheme was used for Greenland simulations in Lipscomb et al. (2019, GMD)
          !  and is retained for back compatibility.
          ! Note: Instead of using a linear ramp for the variation between overburden and the small value
          !       (as for the BPMP and BWAT options above), we use the published formulation of Bueler & van Pelt (2015).
          !       This formulation has N = P_0 for bwat up to ~0.6*effecpress_bwat_threshold; then N decreases
          !        as bwat => effecpress_bwat_threshold.
          !       See Fig. 1b of Bueler & van Pelt (2015).
          ! Note: relative bwat used to be computed in terms of basal_hydro%bwat_till_max.
          !       This formulation gives the same answer, provided that effecpress_bwat_threshold = bwat_till_max.
          !       Both parameters have default values of 2 m.

          do j = 1, nsn
             do i = 1, ewn

                if (bwat(i,j) > 0.0d0) then

                   relative_bwat = max(0.0d0, min(bwat(i,j)/basal_physics%effecpress_bwat_threshold, 1.0d0))

                   ! Eq. 23 from Bueler & van Pelt (2015)
                   basal_physics%effecpress(i,j) = basal_hydro%N_0  &
                        * (basal_physics%effecpress_delta * overburden(i,j) / basal_hydro%N_0)**relative_bwat  &
                        * 10.d0**((basal_hydro%e_0/basal_hydro%C_c) * (1.0d0 - relative_bwat))

                   ! The following line (if uncommented) would implement Eq. 5 of Aschwanden et al. (2016).
                   ! Results are similar to Bueler & van Pelt, but the dropoff in N from P_0 to delta*P_0 begins
                   !  with a larger value of bwat (~0.7*bwat_threshold instead of 0.6*bwat_threshold).

!!                 basal_physics%effecpress(i,j) = basal_physics%effecpress_delta * overburden(i,j)  &
!!                      * 10.d0**((basal_hydro%e_0/basal_hydro%C_c) * (1.0d0 - relative_bwat))

                   ! limit so as not to exceed overburden
                   basal_physics%effecpress(i,j) = min(basal_physics%effecpress(i,j), overburden(i,j))
                end if
             enddo
          enddo

       endif   ! present(bwat)

       !TODO - Not sure this is needed, because of weighting by f_ground.
       where (floating_mask == 1)
          ! set to zero for floating ice
          basal_physics%effecpress = 0.0d0
       end where

    end select   ! which_effecpress

    ! Optionally, reduce N for ice grounded below sea level based on connectivity of subglacial water to the ocean.
    ! N is weighted by the factor (1 - Hf/H)^p, where Hf is the flotation thickness.
    ! p = 1 => full connectivity
    ! 0 < p < 1 => partial connectivity
    ! p = 0 => no connectivity; p_w = 0

    ocean_p = basal_physics%p_ocean_penetration

    effecpress_ocean_p(:,:) = overburden(:,:)

    if (ocean_p > 0.0d0) then

       ! Compute N as a function of f_pattyn = -rhoo*(tops-eus) / (rhoi*thck)
       !   f_pattyn < 0 for land-based ice, < 1 for grounded ice, = 1 at grounding line, > 1 for floating ice
       !TODO - Try averaging thck and topg to vertices, and computing f_pattyn based on these averages?
       !       Might not be as dependent on whether neighbor cells are G or F.

       do j = 1, nsn
          do i = 1, ewn
             if (thck(i,j) > 0.0d0) then
                f_pattyn = rhoo*(eus-topg(i,j)) / (rhoi*thck(i,j))     ! > 1 for floating, < 1 for grounded
                f_pattyn_capped = max( min(f_pattyn, 1.0d0), 0.0d0)    ! capped to lie in the range [0,1]
                effecpress_ocean_p(i,j) = overburden(i,j) * (1.0d0 - f_pattyn_capped)**ocean_p
             endif
          enddo
       enddo

       !WHL - debug
       if (present(itest) .and. present(jtest) .and. present(rtest)) then
          if (this_rank == rtest .and. verbose_effecpress) then

             ! Compute f_pattyn as a 2D field for diagnostics.
             do j = 1, nsn
                do i = 1, ewn
                   if (thck(i,j) > 0.0d0) then
                      f_pattyn_2d(i,j) = rhoo*(eus-topg(i,j)) / (rhoi*thck(i,j))    ! > 1 for floating, < 1 for grounded
                   else  ! no ice
                      if (topg(i,j) - eus >= 0.0d0) then  ! ice-free land
                         f_pattyn_2d(i,j) = 0.0d0
                      else  ! ice-free ocean
                         f_pattyn_2d(i,j) = 1.0d0
                      endif
                   endif
                enddo
             enddo

             print*, ' '
             print*, 'f_pattyn, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.4)',advance='no') f_pattyn_2d(i,j)
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'multiplier for N, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   f_pattyn_capped = max( min(f_pattyn_2d(i,j), 1.0d0), 0.0d0)
                   write(6,'(f10.4)',advance='no') (1.0d0 - f_pattyn_capped)**ocean_p
                enddo
                write(6,*) ' '
             enddo
             print*, ' '
             print*, 'N_ocean_p, itest, jtest, rank =', itest, jtest, rtest
             do j = jtest+3, jtest-3, -1
                write(6,'(i6)',advance='no') j
                do i = itest-3, itest+3
                   write(6,'(f10.0)',advance='no') effecpress_ocean_p(i,j)
                enddo
                write(6,*) ' '
             enddo
          endif   ! verbose_effecpress
       endif   ! present(itest,jtest,rtest)

    else   ! ocean_p = 0

       ! do nothing, (1 - Hf/H)^p = 1

       ! Note(WHL): If ocean_p = 0, then we have N = rhoi*grav*H for floating ice (f_pattyn_capped = 1).
       !            Equivalently, we are defining 0^0 = 1 for purposes of the Leguy et al. effective pressure parameterization.
       !            This is OK for the Schoof basal friction law provided that the resulting beta is multiplied by f_ground,
       !             where f_ground is the fraction of floating ice at a vertex, with f_ground = 0 if all four neighbor cells are floating.
       !            If we were to set N = 0 where f_pattyn_capped = 1 (i.e., defining 0^0 = 0), then we would have a
       !             sudden sharp increase in N_stag (the effective pressure at the vertex) when f_pattyn_capped at a cell center
       !             falls from 1 to a value slightly below 1.  This sudden increase would occur despite the use of a GLP.

    endif

    ! Choose the minimum of the ocean-connection value and the previously computed value.
    ! Thus, the effective pressure can be reduced by an ocean connection or by the presence of meltwater,
    !  but these two processes do not compound on each other.

    basal_physics%effecpress = min(basal_physics%effecpress, effecpress_ocean_p)

    ! Cap the effective pressure at 0x and 1x overburden pressure to avoid strange values going to the friction laws.
    ! This capping may not be necessary, but is included as a precaution.

    where (basal_physics%effecpress < 0.0d0)
       basal_physics%effecpress = 0.0d0
    elsewhere (basal_physics%effecpress > overburden)
       basal_physics%effecpress = overburden
    endwhere

    ! Halo update before staggering
    call parallel_halo(basal_physics%effecpress, parallel)

    ! Interpolate the effective pressure to the staggered grid.
    ! stagger_margin_in = 0: Interpolate using values in all cells, including ice-free cells
    ! (to give a smooth transition in N_stag as a cell switches from ice-free to ice-covered)
    !TODO - Does ice_mask need to be passed in? Modify glissade_stagger so it can be called without a mask.

    call glissade_stagger(ewn,                       nsn,                             &
                          basal_physics%effecpress,  basal_physics%effecpress_stag,   &
                          ice_mask,                  stagger_margin_in = 0)

    if (verbose_effecpress .and. this_rank == rtest) then
       print*, ' '
       print*, 'ocean_p N/overburden, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
             if (thck(i,j) > 0.0d0) then
                write(6,'(f10.5)',advance='no') effecpress_ocean_p(i,j) / overburden(i,j)
             else
                write(6,'(f10.5)',advance='no') 0.0d0
             endif
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'Final N/overburden, itest, jtest, rank =', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
         do i = itest-3, itest+3
             if (overburden(i,j) > 0.0d0) then
                write(6,'(f10.5)',advance='no') basal_physics%effecpress(i,j) / overburden(i,j)
             else
                write(6,'(f10.5)',advance='no') 0.0d0
             endif
          enddo
          write(6,*) ' '
       enddo
       print*, ' '
       print*, 'effecpress_stag:'
       do j = jtest+3, jtest-3, -1
          write(6,'(i6)',advance='no') j
          do i = itest-3, itest+3
            write(6,'(f10.0)',advance='no') basal_physics%effecpress_stag(i,j)
          enddo
          write(6,*) ' '
       enddo
    endif

  end subroutine calc_effective_pressure

!***********************************************************************

  subroutine set_coulomb_c_elevation(ewn,        nsn,   &
                                     topg,       eus,   &
                                     basal_physics,     &
                                     coulomb_c)

    ! Compute coulomb_c as a function of bed elevation.
    ! Assume a linear ramp between the max value at elevation bedmax and the min value at bedmin.

    use glissade_grid_operators, only: glissade_stagger

    integer, intent(in) :: &
         ewn, nsn            ! grid dimensions

    real(dp), dimension(ewn,nsn), intent(in)      :: topg            ! bed topography (m)
    real(dp), intent(in)                          :: eus             ! eustatic sea level (m) relative to z = 0
    type(glide_basal_physics), intent(in)         :: basal_physics   ! basal physics object
    real(dp), dimension(ewn-1,nsn-1), intent(out) :: coulomb_c       ! 2D field of coulomb_c

    real(dp), dimension(ewn-1,nsn-1) :: &
         stagtopg                             ! topg (m) on the staggered grid

    real(dp) :: coulomb_c_min, coulomb_c_max  ! min and max values of coulomb_c (unitless);
                                              ! analogous to tan(phimin) and tan(phimax)
    real(dp) :: bedmin, bedmax                ! bed elevations (m) below which coulomb_c = coulomb_c_min
                                              !  and above which coulomb_c = coulomb_c_max

    real(dp) :: bed                           ! bed elevation (m)
    integer :: ew, ns

    coulomb_c_min = basal_physics%coulomb_c_min
    coulomb_c_max = basal_physics%coulomb_c_max
    bedmin = basal_physics%coulomb_c_bedmin
    bedmax = basal_physics%coulomb_c_bedmax

    ! Interpolate topg to the staggered grid
    ! stagger_margin_in = 0: Interpolate using values in all cells, including ice-free cells

    call glissade_stagger(ewn,         nsn,         &
                          topg,        stagtopg,    &
                          stagger_margin_in = 0)

    ! Compute coulomb_c based on bed elevation
    do ns = 1, nsn-1
       do ew = 1, ewn-1
          bed = stagtopg(ew,ns) - eus
          if (bed <= bedmin) then
             coulomb_c(ew,ns) = coulomb_c_min
          elseif (bed >= bedmax) then
             coulomb_c(ew,ns) = coulomb_c_max
          else   ! bed elevation is between bedmin and bedmax
             coulomb_c(ew,ns) = coulomb_c_min + &
                  ((bed - bedmin)/(bedmax - bedmin)) * (coulomb_c_max - coulomb_c_min)
          endif
       enddo
    enddo

  end subroutine set_coulomb_c_elevation

!=======================================================================

end module glissade_basal_traction

!=======================================================================
