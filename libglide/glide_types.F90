!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_types.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glide_types

  !> Holds type definitions for the derived types used by each 
  !> instance of the ice model. Originally, each of these types
  !> was a module containing variables, which were used as containers
  !> for global variables. However, the need to allow for multiple
  !> ice model instances meant that the nested derived types were instituted
  !> instead. However, there is probably one too many levels in this scheme. 
  !> It would be better if the different types here were contained in the 
  !> higher-level instance type (\texttt{glint\_instance}), rather than 
  !> the intermediate model type (\texttt{glide\_global\_type}). 
  !> 
  !> Note that this \emph{is} now where the defaults are defined for these
  !> variables.

!TODO - Clean up the glide_global type so it holds fewer subtypes?
!       For example, we could replace some work types (tempwk, velowk) with local arrays and parameters.

  use glimmer_global, only: sp, dp, fname_length
  use glimmer_physcon, only: rhoi, rhoo, coni
  use glimmer_paramets, only: unphys_val
  use glimmer_ncdf, only: glimmer_nc_input, glimmer_nc_output, glimmer_funits, &
       NO_RESTART, STANDARD_RESTART, HYBRID_RESTART
  use profile, only: profile_type
  use glimmer_coordinates, only: coordsystem_type
  use glimmer_map_types, only: glimmap_proj
  use glimmer_sparse_type, only: sparse_matrix_type, &
       HO_SPARSE_PCG_INCH, HO_SPARSE_BICG, HO_SPARSE_GMRES, &
       HO_SPARSE_PCG_STANDARD, HO_SPARSE_PCG_CHRONGEAR
  use cism_parallel, only: parallel_type

  implicit none

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Constants that describe the options available
  ! We use these integer parameters elsewhere in the code to avoid
  !  hardwiring of option numbers

  ! basic Glimmer/Glide options

  integer, parameter :: DYCORE_GLIDE = 0        ! old shallow-ice dycore from Glimmer
!WHL - The Glam dycore is no longer supported; choosing whichdycore = 1 will give a fatal error.
!      Might repurpose the '1' slot for a future version of Glissade.
  integer, parameter :: DYCORE_GLAM = 1         ! Payne-Price finite-difference solver

  integer, parameter :: EVOL_PSEUDO_DIFF = 0    ! glide only
  integer, parameter :: EVOL_ADI = 1            ! glide only
  integer, parameter :: EVOL_DIFFUSION = 2      ! glide only
  integer, parameter :: EVOL_INC_REMAP = 3      ! glissade only
  integer, parameter :: EVOL_UPWIND = 4         ! glissade only
  integer, parameter :: EVOL_NO_THICKNESS = 5   ! glissade only

  !NOTE: Use option 1 for prognostic temperature with any dycore
  !      Option 3 is under construction

  integer, parameter :: TEMP_SURFACE_AIR_TEMP = 0
  integer, parameter :: TEMP_PROGNOSTIC = 1
  integer, parameter :: TEMP_STEADY = 2
  integer, parameter :: TEMP_ENTHALPY = 3

  integer, parameter :: TEMP_INIT_ZERO = 0
  integer, parameter :: TEMP_INIT_ARTM = 1
  integer, parameter :: TEMP_INIT_LINEAR = 2

  integer, parameter :: FLWA_CONST_FLWA = 0
  integer, parameter :: FLWA_PATERSON_BUDD_CONST_TEMP = 1
  integer, parameter :: FLWA_PATERSON_BUDD = 2

  integer, parameter :: BTRC_CONSTANT = 1
  integer, parameter :: BTRC_CONSTANT_BWAT = 2
  integer, parameter :: BTRC_CONSTANT_BPMP = 3
  integer, parameter :: BTRC_LINEAR_BMLT = 4
  integer, parameter :: BTRC_TANH_BWAT = 5

  ! Old Glimmer basal water options
  integer, parameter :: BWATER_NONE  = 0
  integer, parameter :: BWATER_LOCAL = 1
  integer, parameter :: BWATER_FLUX  = 2
  integer, parameter :: BWATER_CONST = 3

  integer, parameter :: BMLT_FLOAT_EXTERNAL = 4

  !TODO - Deprecate the quadratic option?

  integer, parameter :: SMB_INPUT_MYR_ICE = 0     ! use 'acab' for input
  integer, parameter :: SMB_INPUT_MMYR_WE = 1     ! use 'smb' for input

 

  integer, parameter :: GTHF_UNIFORM = 0
  integer, parameter :: GTHF_COMPUTE = 2

  integer, parameter :: RELAXED_TOPO_INPUT = 1    ! set relx to input topg
  integer, parameter :: RELAXED_TOPO_COMPUTE = 2  ! Input topo in isostatic equilibrium
                                                  ! compute relaxed topo

  integer, parameter :: ISOSTASY_COMPUTE = 1

  integer, parameter :: LITHOSPHERE_LOCAL = 0
  integer, parameter :: LITHOSPHERE_ELASTIC = 1

  integer, parameter :: ASTHENOSPHERE_FLUID = 0
  integer, parameter :: ASTHENOSPHERE_RELAXING_CONST = 1
  integer, parameter :: ASTHENOSPHERE_RELAXING_LATVAR = 2

  integer, parameter :: CALVING_NONE = 0
  integer, parameter :: CALVING_FLOAT_ZERO = 1
  integer, parameter :: CALVING_FLOAT_FRACTION = 2
  integer, parameter :: CALVING_RELX_THRESHOLD = 3
  integer, parameter :: CALVING_TOPG_THRESHOLD = 4
  integer, parameter :: CALVING_GRID_MASK = 5
  integer, parameter :: CF_ADVANCE_RETREAT_RATE = 6
  integer, parameter :: CALVING_THCK_THRESHOLD = 7
  integer, parameter :: CALVING_STRESS = 8
  integer, parameter :: CALVING_STRESS_STOCHASTIC = 9
  integer, parameter :: CALVING_DAMAGE = 10
  integer, parameter :: EIGEN_CALVING = 11
  integer, parameter :: CALVING_HUYBRECHTS = 12

  integer, parameter :: CALVING_INIT_ON = 1

  integer, parameter :: VERTINT_STANDARD = 0
  integer, parameter :: VERTINT_KINEMATIC_BC = 1

  integer, parameter :: DT_STEPS_PER_YEAR = 1

  integer, parameter :: DM_DT_DIAG_KG_S = 0
  integer, parameter :: DM_DT_DIAG_GT_Y = 1

  integer, parameter :: DIAG_MINTHCK_ZERO = 0
  integer, parameter :: DIAG_MINTHCK_THKLIM = 1

  integer, parameter :: SIGMA_COMPUTE_GLIDE = 0
  integer, parameter :: SIGMA_EXTERNAL = 1
  integer, parameter :: SIGMA_CONFIG = 2
  integer, parameter :: SIGMA_COMPUTE_EVEN = 3
  integer, parameter :: SIGMA_COMPUTE_PATTYN = 4

  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_general

    !> Holds fundamental parameters of the ice model geometry.

    integer :: ewn = 0  !> The number of grid-points in the E-W direction.
    integer :: nsn = 0  !> The number of grid-points in the N-S direction.
    integer :: upn = 3  !> The number of vertical levels in the model.

    type(coordsystem_type) :: ice_grid  !> coordinate system of the ice grid
    type(coordsystem_type) :: velo_grid !> coordinate system of the velocity grid

    ! horizontal dimensions: (x1,y1) for unstaggered ice grid, (x0,y0) for staggered velocity grid
    real(dp), dimension(:),pointer :: x0 => null() !original x0 grid 
    real(dp), dimension(:),pointer :: y0 => null() !original y0 grid
    real(dp), dimension(:),pointer :: x1 => null() !original x1 grid
    real(dp), dimension(:),pointer :: y1 => null() !original y1 grid

    ! global versions of horizontal dimension arrays
    ! These contain all values in the global domain, not just those on the local block
    real(dp), dimension(:),pointer :: x0_global => null()
    real(dp), dimension(:),pointer :: y0_global => null()
    real(dp), dimension(:),pointer :: x1_global => null()
    real(dp), dimension(:),pointer :: y1_global => null()

    ! latitude and longitude on the unstaggered ice grid
    real(dp), dimension(:,:),pointer :: lat => null()
    real(dp), dimension(:,:),pointer :: lon => null()

    ! added to handle the active-blocks option
    integer, dimension(:,:), pointer :: ice_domain_mask => null()  ! = 1 for cells that are potentially active

  end type glide_general

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_options

    !> Holds user options controlling the methods used in the ice-model integration.

    !-----------------------------------------------------------------------
    ! standard options
    !-----------------------------------------------------------------------

    integer :: whichdycore = 0

    ! Choice of two dycores:
    !> \begin{description} 
    !> \item[0] Glide dycore (SIA, serial (SLAP) only)
    !> \item[1] Glam dycore (no longer supported!)
    !> \item[2] Glissade dycore (1st-order, FEM, serial (SLAP) or parallel (F90 native PCG solver) )
    !> \end{description}

    integer :: whichevol = 0

    !TODO: Break into two options: one set (0, 1, 2) for Glide and one set (3, 4, 5) for Glissade
    !> Thickness evolution method:
    !> \begin{description}
    !> \item[0] Pseudo-diffusion 
    !> \item[1] Alternating direction implicit (ADI)
    !> \item[2] Diffusion (also calculates velocities) 
    !> \item[3] Incremental remapping
    !> \item[4] 1st-order upwind scheme
    !> \item[5] Temperature evolves but thickness does not
    !> \end{description}

    integer :: whichtemp = 1

    !> Method of ice temperature calculation:
    !> \begin{description} 
    !> \item[0] Set column to surface air temperature
    !> \item[1] Prognostic temperature solution 
    !> \item[2] Do NOTHING - hold temperatures steady at initial value  
    !> \item[3] Prognostic enthalpy solution
    !> \end{description}

    integer :: temp_init = 1

    ! Temperature initialization:
    !> \begin{description} 
    !> \item[0] Initialize temperature to 0 C
    !> \item[1] Initialize temperature to surface air temperature
    !> \item[2] Initialize temperature with a linear profile in each column
    !> \item[3] Initialize temperature with an advective-diffusive balance in each column
    !> \item[4] Initialize temperature from external file
    !> \end{description}

    !> Method for calculating flow factor $A$:

    integer :: whichflwa = 2

    !> \begin{description} 
    !> \item[0] Set equal to $1\times 10^{-16} \mathrm{Pa}^{-n} \mathrm{yr}^{-1}$
    !> \item[1] \emph{Paterson and Budd} relationship, 
    !> with temperature set to $-5^{\circ}\mathrm{C}$ 
    !> \item[2] \emph{Paterson and Budd} relationship
    !> \end{description}

    integer :: whichbtrc = 0

    !> Basal slip coefficient (used for glide and for glissade local SIA; not glissade HO):
    !> \begin{description}
    !> \item[0] Set equal to zero everywhere
    !> \item[1] Set to (non--zero) constant
    !> \item[2] Set to (non--zero) constant where basal water is present, otherwise to zero
    !> \item[3] Set to (non--zero) constant where bed temperature is at pressure melting point, otherwise to zero
    !> \item[4] linear function of basal melt rate
    !> \item[5] \texttt{tanh} function of basal water depth 
    !> \end{description}

    integer :: whichbwat = 0

    !> Basal water depth: 
    !> \begin{description} 
    !> \item[0] Set to zero everywhere 
    !> \item[1] Compute from local basal water balance 
    !> \item[2] Compute the basal water flux, then find depth via calculation
    !> \item[3] Set to constant (10 m) everywhere, to force T = Tpmp.
    !> \end{description}

    integer :: basal_mbal = 1

    !> basal mass balance:
    !> \begin{description}
    !> \item[0] Basal mass balance not included in continuity equation
    !> \item[1] Basal mass balance included in continuity equation
    !> \end{description}

    !TODO - Change name to smb_input_units?
    integer :: smb_input = 0

    !> units for SMB input:
    !> \begin{description}
    !> \item[0] SMB input in units of m/yr ice (same as acab)
    !> \item[1] SMB input in units of mm/yr water equivalent
    !> \end{description}
    

    integer :: gthf = 0

    !> geothermal heat flux:
    !> \begin{description}
    !> \item[0] prescribed uniform geothermal heat flux
    !> \item[1] read 2D geothermal flux field from input file (if present)
    !> \item[2] calculate geothermal flux using 3d diffusion
    !> \end{description}

    integer :: isostasy = 0

    !> isostasy:
    !> \begin{description}
    !> \item[0] no isostatic adjustment
    !> \item[1] compute isostatic adjustment using lithosphere/asthenosphere model
    !> \end{description}

    integer :: whichcalving = 1

    !> Calving: 
    !> \begin{description} 
    !> \item[0] No calving 
    !> \item[1] Set thickness to zero if floating
    !> \item[2] Lose a fraction of floating ice at marine margin
    !> \item[3] Set thickness to zero if relaxed bedrock is more than a
    !>          certain water depth (variable "marine_limit" in glide_types)  
    !> \item[4] Set thickness to zero if present bedrock topography lies below
    !>          a certain water depth (variable "marine_limit" in glide_types)  
    !> \item[5] Calve based on grid location using a prescribed mask
    !>          (field 'calving_mask' or 'subgrid_calving_mask')
    !> \item[6] Prescribe the rate of calving front advance or retreat
    !> \item[7] Calve ice whose thickness is below a given threshold
    !> \item[8] Deterministic calving based on eigenvalues of the horizontal stress tensor
    !> \item[9] Stochastic calving based on eigenvalues of the horizontal stress tensor
    !> \item[10] Calving based on stochastic internal damage
    !> \item[11] Eigencalving following Levermann et al. (2012)
    !> \item[12] Huybrechts calving
    !TODO - Revise the numbering?
    !> \end{description}

    integer :: calving_init = 0
    !> \begin{description}
    !> \item[0] Do not calve at initialization
    !> \item[1] Calve at initialization
    !> \end{description}

    !Note: This used to be called dt_option; renamed to avoid a naming conflict with CESM
    integer :: dt_input_option = 0
    !> \begin{description}
    !> \item[0] Input dt in years
    !> \item[1] Input dt in steps per year
    !> \end{description}

    integer :: dm_dt_diag = 0
    !> \begin{description}
    !> \item[0] Write dmass/dt diagnostic in units of kg/s
    !> \item[1] Write dmass/dt diagnostic in units of Gt/yr
    !> \end{description}

    integer :: diag_minthck = 1
    !> \begin{description}
    !> \item[0] Include cells with H > 0 in global diagnostics
    !> \item[1] Include cells with H > thklim in global diagnostics
    !> \end{description}

    integer :: whichwvel = 0

    !> Vertical velocities: 
    !> \begin{description}
    !> \item[0] Usual vertical integration 
    !> \item[1] Vertical integration constrained to obey upper kinematic BC
    !> \end{description}

    integer :: which_sigma = 0

    !> \begin{description}
    !> \item[0] compute standard Glimmer sigma coordinates
    !> \item[1] sigma coordinates are given in external file
    !> \item[2] sigma coordinates are given in configuration file
    !> \item[3] evenly spaced levels
    !> \item[4] compute Pattyn sigma coordinates
    !> \end{description}

    logical :: forcewrite_final = .false.
    !> if true, then force a write to output and restart files when the model finishes

    !TODO - Change 'is_restart' to 'restart'
    integer :: is_restart = 0
    !> if the run is a restart of a previous run
    !> \begin{description}
    !> \item[0] normal start-up (take init fields from .nc input file OR if absent, use default options)
    !> \item[1] restart model from previous run (do not calc. temp, rate factor, or vel)
    !> \item[2] hybrid restart; use restart from previous run as the input file, and reset the time
    !> \end{description}

    ! This is a Glimmer serial option
    ! The parallel code enforces periodic EW and NS boundary conditions by default
    logical :: periodic_ew = .false.

    !> \begin{description}
    !> \item[0] no periodic EW boundary conditions
    !> \item[1] periodic EW boundary conditions
    !> \end{description}

  end type glide_options

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_geometry

    !> Holds fields and other information relating to the geometry of the ice sheet and bedrock.

    real(dp),dimension(:,:),pointer :: thck => null()
    !> ice thickness (m)

    real(dp),dimension(:,:),pointer :: usrf => null()
    ! elevation of the upper ice surface (m)

    real(dp),dimension(:,:),pointer :: lsrf => null() 
    ! elevation of the lower ice surface (m)

    real(dp),dimension(:,:),pointer :: topg => null() 
    !> elevation of the bed topography (m)

    real(dp),dimension(:,:,:),pointer :: ice_age => null()
    !> The age of a given ice layer; units are s internally, yr for I/O

    real(dp),dimension(:,:),pointer :: cell_area => null()          !> grid cell area (m^2)

    integer, dimension(:,:),pointer :: thkmask => null()
    !> see glide_mask.f90 for possible values

    ! fields for ice removal, e.g., removal of ice caps
    real(dp),dimension(:,:),pointer :: removal_thck => null()        !> thickness loss of ice removed in grid cell

    integer, dimension(:,:),pointer :: thck_index => null()
    ! Set to nonzero integer for ice-covered cells (thck > 0), cells adjacent to ice-covered cells,
    !  and cells with acab > 0.  The non-zero points are numbered in sequence from the bottom left 
    !  to the top right, going along the rows.

    integer :: totpts = 0       ! total number of points with nonzero thck_index
    logical :: empty = .true.   ! true if totpts = 0

  end type glide_geometry

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_scalars

    !> Holds various global diagnostic scalars that can be written to output

    ! global scalars
    real(dp) :: iarea                  ! total ice area (m^2)
    real(dp) :: iareag                 ! total grounded ice area (m^2)
    real(dp) :: iareaf                 ! total floating ice area (m^2)
    real(dp) :: ivol                   ! total ice volume (m^3)
    real(dp) :: ivol_above_flotation   ! total ice volume above flotation (m^3)
    real(dp) :: imass                  ! total ice mass (kg)
    real(dp) :: imass_above_flotation  ! total ice mass above flotation (kg)

 end type glide_scalars

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_geomderv

    !> Holds the horizontal and temporal derivatives of the thickness and
    !> upper surface elevation, as well as the thickness on the staggered grid.

    !First derivatives on a staggered grid
    real(dp),dimension(:,:),pointer :: dthckdew => null() !> E-W derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdew => null() !> E-W derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdns => null() !> N-S derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdns => null() !> N-S derivative of upper surface elevation.

    !Time derivatives
    real(dp),dimension(:,:),pointer :: dthckdtm => null() !> Temporal derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdtm => null() !> Temporal derivative of upper surface elevation.

    !TODO - Move staggered variables from glide_geomderv type to glide_geometry?

    !Staggered grid versions of geometry variables
    real(dp),dimension(:,:),pointer :: stagthck => null() !> Thickness averaged onto the staggered grid.

  end type glide_geomderv

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_velocity

    !> Holds the velocity fields in 2D and 3D. Some of these fields are defined only in Glide.
    real(dp),dimension(:,:,:),pointer :: uvel  => null()   !> 3D $x$-velocity.
    real(dp),dimension(:,:,:),pointer :: vvel  => null()   !> 3D $y$-velocity.
    real(dp),dimension(:,:,:),pointer :: velnorm => null() ! horizontal ice speed
    real(dp),dimension(:,:,:),pointer :: wvel  => null()   !> 3D $z$-velocity.
    real(dp),dimension(:,:,:),pointer :: wgrd  => null()   !> 3D grid vertical velocity.
    real(dp),dimension(:,:)  ,pointer :: uflx  => null()   !> 
    real(dp),dimension(:,:)  ,pointer :: vflx  => null()   !> 
    real(dp),dimension(:,:)  ,pointer :: ubas  => null()   !> basal $x$-velocity at cell vertices (Glide only)
    real(dp),dimension(:,:)  ,pointer :: vbas  => null()   !> basal $y$-velocity at cell vertices (Glide only)
    real(dp),dimension(:,:)  ,pointer :: diffu => null()   !> 
    real(dp),dimension(:,:)  ,pointer :: total_diffu => null() !> total diffusivity

    ! Note: Velocities have units of m/s within the code.
    ! For I/O, velocities have units of m/yr, with a scale factor of 1/scyr

    real(dp),dimension(:,:)  ,pointer :: btrc  => null()        !> basal traction (scalar field)
    real(dp),dimension(:,:)  ,pointer :: beta  => null()        !> basal shear coefficient on velo grid (Pa yr/m)
    real(dp),dimension(:,:)  ,pointer :: tau_x => null()        !> SIA basal shear stress, x-dir
    real(dp),dimension(:,:)  ,pointer :: tau_y => null()        !> SIA basal shear stress, y-dir

  end type glide_velocity

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Make eus a config file parameter.
  ! Note: smb and acab are both SMB variables, but they have different units.
  !       smb has units of kg/m^2/yr w.e. = mm/yr w.e.
  !       acab has units of m/yr ice for I/O, and SI units in most of the code.
  !       Terms with a prefix 'smb_' or 'acab_' follow the same units convention.

  ! Note on acab_applied_tavg: This is the average value of acab_applied over an output interval.
  !                    If 'average = 1' in the acab_applied entry of glide_vars.def, then acab_applied_tavg is automatically
  !                     accumulated and averaged during runtime, without any additional code needed.
  !                    Other variables with a '_tavg' suffix are handled similarly.
  !                    Between writes, the array acab_applied_tavg holds the running sum of acab_applied*tinc, not the average.
  !                    The sum is divided by the total time only when the output variable acab_applied_tavg is written,
  !                     and then the array is reset to zero (see NAME_avg_accumulate and NAME_avg_reset
  !                     in ncdf_template.F90.in).
  !                    For exact restart in the middle of an averaging interval, the same array is also written
  !                     to restart files, without dividing by the total time, as the output variable acab_applied_tavg_sum.
  !                     On a standard restart, acab_applied_tavg_sum is read back into the array acab_applied_tavg.
  !                    The '_tavg_sum' variables are generated automatically (generate_ncvars.py) and added
  !                     automatically to the restart variable list for each '_tavg' field in an output file;
  !                     there is no separate '_tavg_sum' array, and users should not list them in output files.
  !
  ! Note on smb_corrected: Optionally, smb can be supplemented with an anomaly.
  !                        The background field, smb, does not include the anomaly.
  !                        Write smb_corrected to the output file to see the modified SMB field.
  ! Note on artm_corrected: Optionally, artm can be supplemented with an anomaly.
  !                         The background field, artm, does not include the anomaly.
  !                         Write artm_corrected to the output file to see the modified artm field.

  type glide_climate
     !> Holds fields used to drive the model
     real(dp),dimension(:,:),pointer :: acab            => null() !> Surface mass balance (m/yr ice for I/O, m/s ice in code)
     real(dp),dimension(:,:),pointer :: acab_applied    => null() !> Surface mass balance applied to ice (m/yr ice)
     real(dp),dimension(:,:),pointer :: acab_applied_tavg => null() !> Surface mass balance applied to ice (m/yr ice, time average)
                                                                  !>    = 0 for ice-free cells with acab < 0
     real(dp),dimension(:,:),pointer :: smb             => null() !> Surface mass balance (mm/yr water equivalent)
     real(dp),dimension(:,:),pointer :: artm            => null() !> Annual mean air temperature (degC)
     real(dp),dimension(:,:),pointer :: artm_corrected  => null() !> Annual mean air temperature with anomaly corrections (degC)
     real(dp),dimension(:,:),pointer :: precip          => null() !> precipitation rate (mm/yr w.e.)

     real(dp) :: eus = 0.d0                         !> eustatic sea level

  end type glide_climate

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_calving
     !> holds fields and parameters related to calving
     real(dp),dimension(:,:),  pointer :: calving_thck => null()   !> thickness loss in grid cell due to calving during one time step (m)
                                                                         !> real instead of integer for use with subgrid calving parameterization
                                                                   !> these cells not allowed to fill until upstream neighbors are full
  
     real(dp) :: marine_limit =  -200.d0         !> value of topg/relx at which floating ice calves (m)
                                                 !> (whichcalving = CALVING_RELX_THRESHOLD, CALVING_TOPG_THRESHOLD)
     real(dp) :: calving_fraction = 0.2d0        !> fractional thickness of floating ice that calves
                                                 !> (whichcalving = CALVING_FLOAT_FRACTION)
                                                 !> WHL - previously defined as the fraction of floating ice that does not calve
                                                 !> if calving_timescale = 0, then the full column calves at once
                                                 !> if used, must be set to a nonzero value in the config file
                                                 !> (whichcalving = CALVING_THCK_THRESHOLD, EIGENCALVING, CALVING_DAMAGE)
                                                 !> 0 = no damage, 1 = total damage (whichcalving = CALVING_DAMAGE)
!     real(dp) :: damage_constant1 = 0.0d0        !> damage constant that multiplies tau_eigen1 (yr^-1)
!     real(dp) :: damage_constant2 = 0.0d0        !> damage constant that multiplies tau_eigen2 (yr^-1)

  end type glide_calving

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_temper

    !> Holds fields relating to temperature.

    !Note: temp, flwa and dissip live on the unstaggered vertical grid
    !       at layer interfaces and have vertical dimension (1:upn).
    !
    !      bheatflx, ucondflx, and lcondflx are defined as positive down,
    !       so they will often be < 0.  
    !      However, bfricflx, bhydroflx and dissipcol are defined to be >= 0.
    !
    !      If bheatflx is read from a data file, be careful about the sign!
    !      In input data, the geothermal heat flux may be defined as positive upward,
    !       whereas bheatflx is defined as positive downward.

    real(dp),dimension(:,:,:),pointer :: temp => null()      !> 3D temperature field.
    real(dp),dimension(:,:),  pointer :: bheatflx => null()  !> basal heat flux (W/m^2) (geothermal, positive down)
    real(dp),dimension(:,:,:),pointer :: flwa => null()      !> Glen's flow factor $A$.
    real(dp),dimension(:,:,:),pointer :: dissip => null()    !> interior heat dissipation rate, divided by rhoi*Ci (deg/s)
    real(dp),dimension(:,:),  pointer :: btemp => null()     !> Basal temperature on ice grid; diagnosed from temp(upn)
    ! Note: With a basal temperature GLP, a grid cell can have different basal temps for grounded and floating ice
    real(dp),dimension(:,:),  pointer :: stagbtemp => null() !> Basal temperature on velo grid
    real(dp),dimension(:,:),  pointer :: bpmp => null()      !> Basal pressure melting point temperature
    real(dp),dimension(:,:),  pointer :: stagbpmp => null()  !> Basal pressure melting point temperature on velo grid
                                                             !> = rhoi * Ci * T for cold ice
    real(dp),dimension(:,:),  pointer :: ucondflx => null()  !> conductive heat flux (W/m^2) at upper sfc (positive down)

    !TODO - Remove some of the following from the derived type
    integer  :: niter   = 0   
    logical  :: newtemps = .false. !> new temperatures

  end type glide_temper

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_basal_melt

     !> Holds fields and parameters relating to basal melting

     ! basal melt fields

     real(dp), dimension(:,:), pointer :: &
          bmlt => null(),                         & !> basal melt rate (> 0 for melt, < 0 for freeze-on)
                                                    !> (m/s); bmlt = bmlt_ground in Glide
          bmlt_applied => null(),                 & !> basal melt rate applied to ice (m/s)
                                                    !>    = 0 for ice-free cells with bmlt > 0
          bmlt_applied_tavg => null(),            & !> basal melt rate applied to ice (m/s, time average)
          bmlt_ground => null()                     !> basal melt rate for grounded ice

  end type glide_basal_melt

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_basal_hydro

     !> Holds variables related to basal hydrology
     !> See glide_bwater.F90 for usage details

     ! Note: Ideally, bwat should have MKS units (m), but currently is scaled.
     real(dp),dimension(:,:),  pointer :: bwat => null()        !> Local basal water depth
     real(dp),dimension(:,:),  pointer :: stagbwat => null()    !> Basal water depth on velo grid
     real(dp),dimension(:,:),  pointer :: bwatflx => null()     !> Basal water flux through a cell per unit cell area (m/s)

  end type glide_basal_hydro

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_lithot_type
     !> holds variables for temperature calculations in the lithosphere

     real(dp),dimension(:,:,:),pointer :: temp => null()    !> Three-dimensional temperature field.
     logical, dimension(:,:), pointer :: mask => null()     !> whether the point has been ice covered at some time

     integer :: num_dim = 1                                 !> either 1 or 3 for 1D/3D calculations

     ! The sparse matrix and linearised arrays
     type(sparse_matrix_type) :: fd_coeff, fd_coeff_slap
     real(dp), dimension(:), pointer :: rhs => null()
     real(dp), dimension(:), pointer :: answer => null()
     real(dp), dimension(:), pointer :: supd => null(), diag => null(), subd => null()

     ! work arrays for solver
     real(dp), dimension(:), pointer :: rwork => null()
     integer, dimension(:), pointer :: iwork => null()
     integer mxnelt

     real(dp), dimension(:), pointer :: deltaz => null()    !> array holding grid spacing in z
     real(dp), dimension(:,:), pointer :: zfactors => null()!> array holding factors for finite differences of vertical diffu
     real(dp) :: xfactor,yfactor !> factors for finite differences of horizontal diffu

     real(dp) :: surft = 2.d0          !> surface temperature, used for calculating initial temperature distribution
     real(dp) :: mart  = 2.d0          !> sea floor temperature 
     integer  :: nlayer = 20           !> number of layers in lithosphere
     real(dp) :: rock_base = -5000.d0  !> depth below sea-level at which geothermal heat gradient is applied
     
     integer :: numt = 0        !> number time steps for spinning up GTHF calculations

     real(dp) :: rho_r = 3300.0d0 !> The density of lithosphere (kg m$^{-3}$)
     real(dp) :: shc_r = 1000.0d0 !> specific heat capcity of lithosphere (J kg$^{-1}$ K$^{-1}$)
     real(dp) :: con_r = 3.3d0    !> thermal conductivity of lithosphere (W m$^{-1}$ K$^{-1}$)

     real(dp) :: diffu = 0. !> diffusion coefficient

  end type glide_lithot_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type isos_elastic
     !> Holds data used by isostatic adjustment calculations

     real(dp) :: d = 0.24d25                !> flexural rigidity (N m)
     real(dp) :: lr                         !> radius of relative stiffness
     real(dp) :: a                          !> radius of disk
     real(dp) :: c1,c2,cd3,cd4              !> coefficients
     real(dp), dimension(:,:), pointer :: w => null() !> matrix operator for lithosphere deformation
     integer :: wsize                       !> size of operator (0:rbel_wsize, 0:rbel_wsize), operator is axis symmetric
  end type isos_elastic

  type isostasy_type
     !> contains isostasy configuration

     integer :: lithosphere = 1
     !> method for calculating equilibrium bedrock depression
     !> \begin{description}
     !> \item[0] local lithosphere, equilibrium bedrock depression is found using Archimedes' principle
     !> \item[1] elastic lithosphere, flexural rigidity is taken into account
     !> \end{description}

     integer :: asthenosphere = 1
     !> method for approximating the mantle
     !> \begin{description}
     !> \item[0] fluid mantle, isostatic adjustment happens instantaneously
     !> \item[1] relaxing mantle, exp adjust toward (relx - load), constant timescale
     !> \item[2] relaxing mantle, exp adjust toward (relx - load), laterally varying timescale
     !> \end{description}

    integer :: which_relaxed = 0

    !> relaxed topography:
    !> \begin{description}
    !> \item[0] get current topo (topg) and relaxed topo (relx) from separate input fields
    !> \item[1] first time slice of input topo is relaxed
    !> \item[2] first time slice of input topo is in isostatic equilibrium
    !> \end{description}

     real(dp) :: tau_relax_const = 3000.d0    ! characteristic time constant of relaxing mantle (yr)
     real(dp) :: load_update_interval = 10.d0 ! update interval for elastic load (yr)
     integer :: nlith                         ! update lithosphere every nlith time steps; tinc * nlith = lithosphere_period
     logical :: new_load = .false.            ! set to true if there is a new surface load
     type(isos_elastic) :: rbel               ! structure holding elastic lithosphere setup

     real(dp),dimension(:,:),pointer :: relx => null()  ! elevation of relaxed topography (m)
     real(dp),dimension(:,:),pointer :: load => null()  ! deflection due to applied load on lithosphere (m)
                                                        ! defined as positive for downward deflection
     real(dp),dimension(:,:),pointer :: load_factors => null() ! temporary used for load calculation

  end type isostasy_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_numerics

    !> Parameters relating to the model numerics
    real(dp) :: tstart =     0.d0 !> starting time
                                  !> Note: For restarts, tstart is set to the time the restart begins
    real(dp) :: tend   =  1000.d0 !> end time
    real(dp) :: time   =     0.d0 !> main loop counter in years
    real(dp) :: tinc   =     1.d0 !> time step of main loop in years 
    real(dp) :: ntem   =     1.d0 !> multiplier of main time step; allows longer temperature time step
    real(dp) :: alpha  =    0.5d0 !> richard suggests 1.5 - was a parameter in original
    real(dp) :: thklim =   100.d0 !> min thickness for computing ice dynamics (m) 
    real(dp) :: dew    =    20.d3     !> grid cell size in east-west direction
    real(dp) :: dns    =    20.d3     !> grid cell size in north-south direction
    real(dp) :: dt     =     0.d0     !> ice dynamics timestep
    real(dp) :: dttem  =     0.d0     !> temperature timestep
    integer  :: nsteps_per_year = 0   !> number of time steps per year (for dt_option = STEPS_PER_YEAR)
    integer  :: tstep_count = 0   !> number of time steps since the start of the simulation
    
    ! Vertical coordinate ---------------------------------------------------
                                                               
    real(dp),dimension(:),pointer :: sigma => null() !> Sigma values for vertical spacing of 
                                                     !> model levels
    real(dp),dimension(:),pointer :: stagsigma => null()     !> Staggered values of sigma (layer midpts)
    real(dp),dimension(:),pointer :: stagwbndsigma => null() !> Staggered values of sigma (layer midpts) with boundaries
    real(dp), dimension(:,:), pointer :: dups => null()      !> vertical grid quantities related to sigma and stagsigma

    integer :: profile_period = 100            ! profile frequency

    real(dp) :: dt_diag = 0.d0     ! diagnostic interval (write diagnostic output every dt_diag years)
                                   ! dt_diag = 0 => never write diagnostic output
    integer  :: ndiag = 0          ! diagnostic interval (write diagnostic output every ndiag timesteps)
                                   ! ndiag = 0 => never write diagnostic output
    integer  :: idiag = 1          ! global grid indices for diagnostic point
    integer  :: jdiag = 1          ! 
    integer  :: idiag_local = 1    ! local grid indices for diagnostic point
    integer  :: jdiag_local = 1
    integer  :: rdiag_local = 0    ! task number for diagnostic point

    ! CFL diagnostics and adaptive subcycling
    real(dp) :: adv_cfl_dt = 0.0d0  ! maximum allowable dt (yrs) based on advective CFL (calculated by model for each time step)
    real(dp) :: diff_cfl_dt = 0.0d0 ! maximum allowable dt (yrs) based on diffusive CFL (calculated by model for each time step)

  end type glide_numerics

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_velowk
    real(dp),dimension(:),  pointer :: depth    => null()
    real(dp),dimension(:),  pointer :: dupsw    => null()
    real(dp),dimension(:),  pointer :: depthw   => null()
    real(dp),dimension(:),  pointer :: suvel    => null()
    real(dp),dimension(:),  pointer :: svvel    => null()
    real(dp),dimension(:,:),pointer :: fslip    => null()
    real(dp),dimension(:,:),pointer :: dintflwa => null()
    real(dp),dimension(:),  pointer :: dups     => null()
    real(dp),dimension(4) :: fact
    real(dp),dimension(4) :: c = 0.d0
    real(dp) :: watwd  = 3.0d0
    real(dp) :: watct  = 10.0d0
!!    real(dp) :: trc0   = 0.d0
    real(dp) :: trcmin = 0.0d0
    real(dp) :: marine = 1.0d0
    real(dp) :: trcmax = 10.0d0
  end type glide_velowk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_thckwk
     real(dp),dimension(:,:),  pointer :: oldthck   => null()
     real(dp),dimension(:,:),  pointer :: oldthck2  => null()
     real(dp),dimension(:,:,:),pointer :: olds      => null()
     integer  :: nwhich  = 2
     real(dp) :: oldtime = 0.d0
     
     ! next four are for ADI evolution only
     real(dp), dimension(:), pointer :: alpha => null()
     real(dp), dimension(:), pointer :: beta  => null()
     real(dp), dimension(:), pointer :: gamma => null()
     real(dp), dimension(:), pointer :: delta => null()

  end type glide_thckwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_tempwk
    real(dp),dimension(:,:,:),pointer :: inittemp => null()
    real(dp),dimension(:,:,:),pointer :: compheat => null()
    real(dp),dimension(:,:,:),pointer :: initadvt => null()
    real(dp),dimension(:),    pointer :: dupa     => null()
    real(dp),dimension(:),    pointer :: dupb     => null()
    real(dp),dimension(:),    pointer :: dupc     => null()
    real(dp),dimension(:),    pointer :: c1       => null()
    real(dp),dimension(:,:),  pointer :: dups     => null()
    real(dp),dimension(:,:),  pointer :: wphi     => null()
    real(dp),dimension(:,:),  pointer :: smth     => null()
    real(dp),dimension(:,:,:),pointer :: hadv_u   => null()
    real(dp),dimension(:,:,:),pointer :: hadv_v   => null()

    !*sfp** added space to the next 2 (cons, f) for use w/ HO and SSA dissip. calc. 
    real(dp),dimension(5)             :: cons     = 0.d0
    real(dp),dimension(5)             :: f        = 0.d0
    real(dp),dimension(8)             :: c        = 0.d0
    real(dp),dimension(2)             :: slide_f
    real(dp),dimension(2) :: advconst = 0.d0
    real(dp) :: zbed        = 0.d0
    real(dp) :: dupn        = 0.d0
    real(dp) :: wmax        = 0.d0
    real(dp) :: dt_wat      = 0.d0
    real(dp) :: watvel      = 0.d0
    integer  :: nwat        = 0
  end type glide_tempwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Move these parameters to types associated with a certain kind of physics
  !TODO - Change default geot to 0.0, so that idealized tests by default have no mass loss
  type glide_paramets
    real(dp),dimension(5) :: bpar = (/ 0.2d0, 0.5d0, 0.0d0 ,1.0d-2, 1.0d0/)
    real(dp) :: btrac_const = 0.d0     ! m yr^{-1} Pa^{-1} (gets scaled during init)
    real(dp) :: btrac_slope = 0.d0     ! Pa^{-1} (gets scaled during init)
    real(dp) :: btrac_max = 0.d0       ! m yr^{-1} Pa^{-1} (gets scaled during init)
    real(dp) :: geot   = -5.0d-2       ! W m^{-2}, positive down  !TODO - Change default to 0.0?
    real(dp) :: flow_enhancement_factor_ground = 1.0d0   ! flow enhancement parameter for the Arrhenius relationship;
                                                         ! grounded ice only; typically > 1 for SIA models to speed up the ice
                                                         ! Default is 1.0, but for marine simulations a smaller value
                                                         !  is often needed to match observed shelf speeds
    real(dp) :: hydtim = 1000.0d0      ! years, converted to s^{-1} and scaled
                                       ! 0 if no drainage
    real(dp) :: bwat_smooth = 0.01d0   ! basal water field smoothing strength

    real(dp) :: default_flwa = 1.0d-16 ! Glen's A to use in isothermal case, in units Pa^{-n} yr^{-1} 
                                       ! (would change to e.g. 4.6e-18 in EISMINT-ROSS case)

  end type glide_paramets

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_prof_type
     integer :: geomderv
     integer :: hvelos
     integer :: ice_mask1
     integer :: temperature
     integer :: ice_evo
     integer :: ice_mask2
     integer :: isos_water
     integer :: isos
  end type glide_prof_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! for the Glide thickness solver
  type, public :: glide_solver

    type(sparse_matrix_type)  :: matrix
    real(dp),dimension(:),pointer :: rhsd    => null()
    real(dp),dimension(:),pointer :: answ    => null()

    !TODO - KJE - Remove ewn/nsn from glide_solver type once new glide_global_type is working and we can use those ewn/nsn
    integer :: ewn
    integer :: nsn

  end type glide_solver

       
  type glide_global_type    ! type containing all of the above for an ice sheet model instance
    integer              :: model_id = 0 !> identifier of this model instance (1..N, where N is the number of ice sheets in this run)
    type(glide_general)  :: general
    type(glide_options)  :: options
    type(glide_geometry) :: geometry
    type(glide_geomderv) :: geomderv
    type(glide_velocity) :: velocity
    type(glide_scalars)  :: scalars
    type(glide_climate)  :: climate
    type(glide_calving)  :: calving
    type(glide_temper)   :: temper
    type(glide_basal_hydro)  :: basal_hydro
    type(glide_basal_melt)   :: basal_melt
    type(glide_lithot_type)  :: lithot
    type(glimmer_funits) :: funits
    type(glide_numerics) :: numerics
    type(glide_velowk)   :: velowk
    type(glide_thckwk)   :: thckwk
    type(glide_tempwk)   :: tempwk
    type(glide_paramets) :: paramets
    type(glimmap_proj)   :: projection
    type(profile_type)   :: profile
    type(glide_prof_type):: glide_prof
    type(isostasy_type)  :: isostasy
    type(glide_solver):: solver_data
    type(parallel_type)  :: parallel  !WHL - new derived type to hold information about parallel communication

  end type glide_global_type

contains

  subroutine glide_allocarr(model)    

    !> Allocates the model arrays, and initialises some of them to zero.
    !> These are the arrays allocated, and their dimensions:

    !TODO - Make sure the itemized lists in subroutine glide_allocarr are complete.

    !> In \texttt{model\%temper}:
    !> \begin{itemize}
    !> \item \texttt{temp(upn,0:ewn+1,0:nsn+1)}   !WHL - 2 choices
    !> \item \texttt{bheatflx(ewn,nsn)}
    !> \item \texttt{flwa(upn,ewn,nsn)}           !WHL - 2 choices
    !> \item \texttt{dissip(upn,ewn,nsn)}         !WHL - 2 choices
    !> \item \texttt{ucondflx(ewn,nsn)}
    !> \end{itemize}

    !> In \texttt{model\%basal_melt}:
    !> \begin{itemize}
    !> \item \texttt{bmlt(ewn,nsn)}
    !> \item \texttt{bmlt_ground(ewn,nsn)}
    !> \item \texttt{bmlt_applied(ewn,nsn)}
    !> \end{itemize}

    !> In \texttt{model\%velocity}:
    !> \begin{itemize}
    !> \item \texttt{uvel(upn,ewn-1,nsn-1))}
    !> \item \texttt{vvel(upn,ewn-1,nsn-1))}
    !> \item \texttt{wvel(upn,ewn,nsn))}
    !> \item \texttt{wgrd(upn,ewn,nsn))}
    !> \item \texttt{uflx(ewn-1,nsn-1))}
    !> \item \texttt{vflx(ewn-1,nsn-1))}
    !> \item \texttt{diffu(ewn,nsn))}
    !> \item \texttt{btrc(ewn-1,nsn-1))}
    !> \end{itemize}

    !> In \texttt{model\%climate}:
    !> \begin{itemize}
    !> \item \texttt{acab(ewn,nsn))}
    !> \item \texttt{artm(ewn,nsn))}
    !> \end{itemize}

    !> In \texttt{model\%geomderv}:
    !> \begin{itemize}
    !> \item \texttt{dthckdew(ewn,nsn))}
    !> \item \texttt{dusrfdew(ewn,nsn))}
    !> \item \texttt{dthckdns(ewn,nsn))}
    !> \item \texttt{dusrfdns(ewn,nsn))}
    !> \item \texttt{dthckdtm(ewn,nsn))}
    !> \item \texttt{dusrfdtm(ewn,nsn))}
    !> \item \texttt{stagthck(ewn-1,nsn-1))}
    !> \end{itemize}
  
    !> In \texttt{model\%geometry}:
    !> \begin{itemize}
    !> \item \texttt{thck(ewn,nsn))}
    !> \item \texttt{usrf(ewn,nsn))}
    !> \item \texttt{lsrf(ewn,nsn))}
    !> \item \texttt{topg(ewn,nsn))}
    !> \item \texttt{mask(ewn,nsn))}
    !* (DFM) added floating_mask, ice_mask, lower_cell_loc, and lower_cell_temp
    !> \end{itemize}

    !> In \texttt{model\%thckwk}:
    !> \begin{itemize}
    !> \item \texttt{olds(ewn,nsn,thckwk\%nwhich))}
    !> \end{itemize}

    !> In \texttt{model\%numerics}:
    !> \begin{itemize}
    !> \item \texttt{sigma(upn))}
    !> \end{itemize}

    !> In \texttt{model\%numerics}:
    !> \begin{itemize}
    !> \item \texttt{stagsigma(upn-1))}
    !> \end{itemize}

    use glimmer_log
    use glimmer_coordinates, only: coordsystem_allocate
    use glimmer_sparse_type, only: new_sparse_matrix

    implicit none

    type(glide_global_type),intent(inout) :: model

    integer :: ewn,nsn,upn               !> local array dimensions
    integer :: global_ewn, global_nsn    !> global array dimensions
    integer :: n

    ! for simplicity, copy these values...
    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn
    global_ewn = model%parallel%global_ewn
    global_nsn = model%parallel%global_nsn
    
    ! horizontal coordinates

    allocate(model%general%x0(ewn-1)); model%general%x0 = 0.d0  ! velocity grid
    allocate(model%general%y0(nsn-1)); model%general%y0 = 0.d0
    allocate(model%general%x1(ewn)); model%general%x1 = 0.d0    ! ice grid (for scalars)
    allocate(model%general%y1(nsn)); model%general%y1 = 0.d0

    allocate(model%general%x0_global(global_ewn-1)); model%general%x0_global = 0.d0  ! velocity grid
    allocate(model%general%y0_global(global_nsn-1)); model%general%y0_global = 0.d0
    allocate(model%general%x1_global(global_ewn)); model%general%x1_global = 0.d0    ! ice grid (for scalars)
    allocate(model%general%y1_global(global_nsn)); model%general%y1_global = 0.d0

    ! vertical sigma coordinates
    ! If we already have sigma, don't reallocate

    if (associated(model%numerics%sigma)) then
       if (size(model%numerics%sigma) /= upn) then
          call write_log('Wrong number of sigma levels given',GM_FATAL)
       end if
    else
       allocate(model%numerics%sigma(upn))
    endif

    allocate(model%numerics%stagsigma(upn-1))
    allocate(model%numerics%stagwbndsigma(0:upn))  !MJH added (0:upn) as separate variable
    allocate(model%numerics%dups(upn+1,2))    !TODO - upn-1 instead?

    ! latitude and longitude
    call coordsystem_allocate(model%general%ice_grid, model%general%lat)
    call coordsystem_allocate(model%general%ice_grid, model%general%lon)

    ! ice domain mask (to identify active blocks)
    call coordsystem_allocate(model%general%ice_grid, model%general%ice_domain_mask)

    !NOTE: The temperature and flow factor live on the unstaggered vertical grid,
    !       and extra rows and columns (with indices 0:ewn+1, 0:nsn+1) are needed.

    allocate(model%temper%temp(upn,0:ewn+1,0:nsn+1))
    call coordsystem_allocate(model%general%ice_grid, upn, model%temper%flwa)
    call coordsystem_allocate(model%general%ice_grid, upn, model%temper%dissip)

    ! MJH - Set temp and flwa to physically unrealistic values so we can tell later if 
    !       arrays were initialized correctly
    model%temper%temp(:,:,:) = unphys_val  ! large negative number
    model%temper%flwa(:,:,:) = unphys_val
    model%temper%dissip(:,:,:) = 0.d0

    call coordsystem_allocate(model%general%ice_grid,  model%temper%bheatflx)
    call coordsystem_allocate(model%general%ice_grid,  model%temper%bpmp)
    call coordsystem_allocate(model%general%velo_grid, model%temper%stagbpmp)
    call coordsystem_allocate(model%general%ice_grid,  model%temper%btemp)
    call coordsystem_allocate(model%general%velo_grid, model%temper%stagbtemp)
    call coordsystem_allocate(model%general%ice_grid,  model%temper%ucondflx)

    call coordsystem_allocate(model%general%ice_grid,  model%basal_hydro%bwat)
    call coordsystem_allocate(model%general%velo_grid, model%basal_hydro%stagbwat)
    call coordsystem_allocate(model%general%ice_grid,  model%basal_hydro%bwatflx)

    ! velocity arrays

    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%uvel)
    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%vvel)
    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%velnorm)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%uflx)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vflx)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%btrc)
    call coordsystem_allocate(model%general%ice_grid,  upn, model%velocity%wvel)

    ! The following are on the extended staggered grid, which is the same size as the ice grid.

    call coordsystem_allocate(model%general%ice_grid,  upn, model%velocity%wgrd)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%ubas)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vbas)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%diffu)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%total_diffu)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%tau_x)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%tau_y)

    ! geometry arrays
    call coordsystem_allocate(model%general%ice_grid, model%geometry%thck)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%usrf)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%lsrf)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%topg)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%thkmask)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%cell_area)

    call coordsystem_allocate(model%general%velo_grid, model%geomderv%stagthck)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dthckdew)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dthckdns)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dusrfdew)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dusrfdns)

    call coordsystem_allocate(model%general%ice_grid, model%geometry%removal_thck)

    ! mass flux arrays

    call coordsystem_allocate(model%general%ice_grid, model%geometry%thck_index)
    call coordsystem_allocate(model%general%ice_grid,  model%geomderv%dthckdtm)
    call coordsystem_allocate(model%general%ice_grid,  model%geomderv%dusrfdtm)
    allocate(model%thckwk%olds(ewn,nsn,model%thckwk%nwhich))
    model%thckwk%olds = 0.0d0
    call coordsystem_allocate(model%general%ice_grid, model%thckwk%oldthck)
    call coordsystem_allocate(model%general%ice_grid, model%thckwk%oldthck2)

    ! bmlt arrays
    call coordsystem_allocate(model%general%ice_grid,  model%basal_melt%bmlt)
    call coordsystem_allocate(model%general%ice_grid,  model%basal_melt%bmlt_applied)
    call coordsystem_allocate(model%general%ice_grid,  model%basal_melt%bmlt_applied_tavg)
    call coordsystem_allocate(model%general%ice_grid,  model%basal_melt%bmlt_ground)

    ! climate arrays
    call coordsystem_allocate(model%general%ice_grid, model%climate%acab)
    call coordsystem_allocate(model%general%ice_grid, model%climate%acab_applied)
    call coordsystem_allocate(model%general%ice_grid, model%climate%acab_applied_tavg)
    call coordsystem_allocate(model%general%ice_grid, model%climate%smb)
    call coordsystem_allocate(model%general%ice_grid, model%climate%artm)
    call coordsystem_allocate(model%general%ice_grid, model%climate%artm_corrected)
    call coordsystem_allocate(model%general%ice_grid, model%climate%precip)

    ! calving arrays
    call coordsystem_allocate(model%general%ice_grid, model%calving%calving_thck)
    if (model%options%whichcalving == CALVING_DAMAGE) then
    else
       ! allocate with size 1, since they need to be allocated to be passed to calving subroutine
    endif

    ! matrix solver arrays
    allocate (model%solver_data%rhsd(ewn*nsn))
    allocate (model%solver_data%answ(ewn*nsn))

    call new_sparse_matrix(ewn*nsn, 5*ewn*nsn, model%solver_data%matrix)

    ! lithosphere arrays

    if (model%options%gthf == GTHF_COMPUTE) then
       allocate(model%lithot%temp(1:ewn,1:nsn,model%lithot%nlayer)); model%lithot%temp = 0.d0
       call coordsystem_allocate(model%general%ice_grid, model%lithot%mask)
    endif

    ! Note: relx needs to be allocated always
    call coordsystem_allocate(model%general%ice_grid, model%isostasy%relx)
    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call coordsystem_allocate(model%general%ice_grid, model%isostasy%load)
       call coordsystem_allocate(model%general%ice_grid, model%isostasy%load_factors)
    endif

  end subroutine glide_allocarr

  subroutine glide_deallocarr(model)

    !> deallocate model arrays
    !TODO - Check that all arrays allocated above are deallocated here.

    implicit none
    type(glide_global_type),intent(inout) :: model

    ! horizontal coordinates

    if (associated(model%general%x0)) &
        deallocate(model%general%x0) 
    if (associated(model%general%y0)) &
        deallocate(model%general%y0) 
    if (associated(model%general%x1)) &
        deallocate(model%general%x1) 
    if (associated(model%general%y1)) &
        deallocate(model%general%y1) 
    if (associated(model%general%x0_global)) &
        deallocate(model%general%x0_global)
    if (associated(model%general%y0_global)) &
        deallocate(model%general%y0_global)
    if (associated(model%general%x1_global)) &
        deallocate(model%general%x1_global)
    if (associated(model%general%y1_global)) &
        deallocate(model%general%y1_global)
    if (associated(model%general%lat)) &
        deallocate(model%general%lat)
    if (associated(model%general%lon)) &
        deallocate(model%general%lon)
    if (associated(model%general%ice_domain_mask)) &
        deallocate(model%general%ice_domain_mask)

    ! vertical sigma coordinates

    if (associated(model%numerics%sigma)) &
        deallocate(model%numerics%sigma)
    if (associated(model%numerics%stagsigma)) &
        deallocate(model%numerics%stagsigma)
    if (associated(model%numerics%stagwbndsigma)) &
        deallocate(model%numerics%stagwbndsigma)

    ! temperature arrays

    if (associated(model%temper%temp)) &
        deallocate(model%temper%temp)
    if (associated(model%temper%bheatflx)) &
        deallocate(model%temper%bheatflx)
    if (associated(model%temper%bpmp)) &
        deallocate(model%temper%bpmp)
    if (associated(model%temper%stagbpmp)) &
        deallocate(model%temper%stagbpmp)
    if (associated(model%temper%btemp)) &
        deallocate(model%temper%btemp)
    if (associated(model%temper%stagbtemp)) &
        deallocate(model%temper%stagbtemp)
    if (associated(model%temper%ucondflx)) &
        deallocate(model%temper%ucondflx)
    if (associated(model%temper%flwa)) &
        deallocate(model%temper%flwa)
    if (associated(model%temper%dissip)) &
        deallocate(model%temper%dissip)

    ! velocity arrays

    if (associated(model%velocity%uvel)) &
        deallocate(model%velocity%uvel)
    if (associated(model%velocity%vvel)) &
        deallocate(model%velocity%vvel)
    if (associated(model%velocity%velnorm)) &
        deallocate(model%velocity%velnorm)
    if (associated(model%velocity%wvel)) &
        deallocate(model%velocity%wvel)
    if (associated(model%velocity%uflx)) &
        deallocate(model%velocity%uflx)
    if (associated(model%velocity%vflx)) &
        deallocate(model%velocity%vflx)
    if (associated(model%velocity%btrc)) &
        deallocate(model%velocity%btrc)
    if (associated(model%velocity%wgrd)) &
        deallocate(model%velocity%wgrd)
    if (associated(model%velocity%ubas)) &
        deallocate(model%velocity%ubas)
    if (associated(model%velocity%vbas)) &
        deallocate(model%velocity%vbas)
    if (associated(model%velocity%diffu)) &
        deallocate(model%velocity%diffu)
    if (associated(model%velocity%total_diffu)) &
        deallocate(model%velocity%total_diffu)
    if (associated(model%velocity%tau_x)) &
        deallocate(model%velocity%tau_x)
    if (associated(model%velocity%tau_y)) &
        deallocate(model%velocity%tau_y)

!!    if (associated(model%velocity%velmask)) &   ! no longer used
!!       deallocate(model%velocity%velmask) 

    if (associated(model%velocity%beta)) &
        deallocate(model%velocity%beta)

    ! basal hydrology arrays
    if (associated(model%basal_hydro%bwat)) &
        deallocate(model%basal_hydro%bwat)
    if (associated(model%basal_hydro%stagbwat)) &
        deallocate(model%basal_hydro%stagbwat)
    if (associated(model%basal_hydro%bwatflx)) &
        deallocate(model%basal_hydro%bwatflx)

    ! basal melt arrays

    if (associated(model%basal_melt%bmlt)) &
        deallocate(model%basal_melt%bmlt)
    if (associated(model%basal_melt%bmlt_applied)) &
        deallocate(model%basal_melt%bmlt_applied)
    if (associated(model%basal_melt%bmlt_applied_tavg)) &
        deallocate(model%basal_melt%bmlt_applied_tavg)
    if (associated(model%basal_melt%bmlt_ground)) &
        deallocate(model%basal_melt%bmlt_ground)

    ! geometry arrays

    if (associated(model%geometry%thck)) &
        deallocate(model%geometry%thck)
    if (associated(model%geometry%usrf)) &
        deallocate(model%geometry%usrf)
    if (associated(model%geometry%lsrf)) &
        deallocate(model%geometry%lsrf)
    if (associated(model%geometry%topg)) &
        deallocate(model%geometry%topg)
    if (associated(model%geometry%thkmask)) &
        deallocate(model%geometry%thkmask)
    if (associated(model%geomderv%stagthck)) &
        deallocate(model%geomderv%stagthck)
    if (associated(model%geomderv%dthckdew)) &
        deallocate(model%geomderv%dthckdew)
    if (associated(model%geomderv%dthckdns)) &
        deallocate(model%geomderv%dthckdns)
    if (associated(model%geomderv%dusrfdew)) &
        deallocate(model%geomderv%dusrfdew)
    if (associated(model%geomderv%dusrfdns)) &
        deallocate(model%geomderv%dusrfdns)
    if (associated(model%geometry%cell_area)) &
        deallocate(model%geometry%cell_area)

    if (associated(model%geometry%removal_thck)) &
       deallocate(model%geometry%removal_thck)

    if (associated(model%geometry%thck_index)) &
        deallocate(model%geometry%thck_index)
    if (associated(model%geomderv%dthckdtm)) &
        deallocate(model%geomderv%dthckdtm)
    if (associated(model%geomderv%dusrfdtm)) &
        deallocate(model%geomderv%dusrfdtm)
    if (associated(model%thckwk%olds)) &
        deallocate(model%thckwk%olds)
    if (associated(model%thckwk%oldthck)) &
        deallocate(model%thckwk%oldthck)
    if (associated(model%thckwk%oldthck2)) &
        deallocate(model%thckwk%oldthck2)

    if (associated(model%geometry%ice_age)) &
        deallocate(model%geometry%ice_age)

    ! climate arrays

    if (associated(model%climate%acab)) &
        deallocate(model%climate%acab)
    if (associated(model%climate%acab_applied)) &
        deallocate(model%climate%acab_applied)
    if (associated(model%climate%acab_applied_tavg)) &
        deallocate(model%climate%acab_applied_tavg)
    if (associated(model%climate%smb)) &
        deallocate(model%climate%smb)
    if (associated(model%climate%artm)) &
        deallocate(model%climate%artm)
    if (associated(model%climate%artm_corrected)) &
        deallocate(model%climate%artm_corrected)
    if (associated(model%climate%precip)) &
        deallocate(model%climate%precip)

    ! calving arrays
    if (associated(model%calving%calving_thck)) &
        deallocate(model%calving%calving_thck)

    ! matrix solver arrays

    if (associated(model%solver_data%rhsd))  &  
        deallocate(model%solver_data%rhsd)
    if (associated(model%solver_data%answ))  &
        deallocate(model%solver_data%answ)

    ! lithosphere arrays

    if (associated(model%lithot%temp)) &
        deallocate(model%lithot%temp)
    if (associated(model%lithot%mask)) &
        deallocate(model%lithot%mask)

    ! isostasy arrays

    if (associated(model%isostasy%relx)) &
        deallocate(model%isostasy%relx)
    if (associated(model%isostasy%load)) &
        deallocate(model%isostasy%load)
    if (associated(model%isostasy%load_factors)) &
        deallocate(model%isostasy%load_factors)

    ! projection arrays
    if (associated(model%projection%stere)) then
       if (associated(model%projection%stere%area_factor)) then
          deallocate(model%projection%stere%area_factor)
       endif
    endif

  end subroutine glide_deallocarr

  ! some accessor functions
  function get_dew(model)
    !> return scaled x node spacing
    implicit none
    real(dp) :: get_dew
    type(glide_global_type) :: model

    get_dew = model%numerics%dew
  end function get_dew

  function get_dns(model)
    !> return scaled y node spacing
    implicit none
    real(dp) :: get_dns
    type(glide_global_type) :: model

    get_dns = model%numerics%dns
  end function get_dns

  function get_tstart(model)
    !> return start time
    implicit none
    real(dp) :: get_tstart
    type(glide_global_type) :: model
    
    get_tstart = model%numerics%tstart
  end function get_tstart

  function get_tend(model)
    !> return end time
    implicit none
    real(dp) :: get_tend
    type(glide_global_type) :: model
    
    get_tend = model%numerics%tend
  end function get_tend

  function get_tinc(model)
    !> return time increment
    implicit none
    real(dp) :: get_tinc
    type(glide_global_type) :: model
    
    get_tinc = model%numerics%tinc
  end function get_tinc

  function get_ewn(model)
    !> get number of nodes in x dir
    implicit none
    integer get_ewn
    type(glide_global_type) :: model

    get_ewn = model%general%ewn
  end function get_ewn

  function get_nsn(model)
    !> get number of nodes in y dir
    implicit none
    integer get_nsn
    type(glide_global_type) :: model

    get_nsn = model%general%nsn
  end function get_nsn

  subroutine set_time(model,time)
    !> Set the model time counter --- useful for
    !> fractional year output
    implicit none
    type(glide_global_type) :: model
    real(dp) :: time

    model%numerics%time = time
  end subroutine set_time

end module glide_types
