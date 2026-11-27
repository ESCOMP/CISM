!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_types.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glissade_types

  !> Holds type definitions for the derived types used by each CISM instance.

  use glimmer_global, only: sp, dp, fname_length
  use glimmer_physcon, only: rhoi, rhoo, coni
  use glimmer_paramets, only: unphys_val
  use glimmer_ncdf, only: glimmer_nc_input, glimmer_nc_output, glimmer_funits, &
       NO_RESTART, STANDARD_RESTART, HYBRID_RESTART
  use cism_profile, only: profile_type
  use glimmer_coordinates, only: coordsystem_type
  use glimmer_map_types, only: glimmap_proj
  use cism_parallel, only: parallel_type
  use glimmer_config, only: option_choice

  implicit none

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Constants that describe the options available
  ! We use these integer parameters elsewhere in the code to avoid hardwiring of option numbers.
  !
  ! Options are given by name in the config file. For each such option, a parameter array of
  ! type option_choice (glimmer_config) lists the integer value, the name and the description of
  ! each allowed value. GetOption reads the option by name, and log_option writes the name and
  ! the description to the log. The constants are named <KEY>_<NAME>.

  ! [grid] global_bc
  integer, parameter :: GLOBAL_BC_PERIODIC = 0  ! doubly periodic
  integer, parameter :: GLOBAL_BC_OUTFLOW = 1   ! free outflow; scalars in global halo set to zero
  integer, parameter :: GLOBAL_BC_NO_PENETRATION = 2   ! no penetration; outflow set to zero at global boundaries
  integer, parameter :: GLOBAL_BC_NO_ICE = 3    ! no ice at global boundary; scalars set to zero adjacent to boundary
  type(option_choice), dimension(4), parameter :: global_bc_choices = [ &
       option_choice(GLOBAL_BC_PERIODIC,       'periodic',       'periodic boundaries'), &
       option_choice(GLOBAL_BC_OUTFLOW,        'outflow',        'outflow at boundaries'), &
       option_choice(GLOBAL_BC_NO_PENETRATION, 'no_penetration', 'no penetration at boundaries'), &
       option_choice(GLOBAL_BC_NO_ICE,         'no_ice',         'no ice at boundaries') ]

  ! [grid] sigma_method
  integer, parameter :: SIGMA_METHOD_GLIDE = 0
  integer, parameter :: SIGMA_METHOD_FILE = 1
  integer, parameter :: SIGMA_METHOD_CONFIG = 2
  integer, parameter :: SIGMA_METHOD_UNIFORM = 3
  integer, parameter :: SIGMA_METHOD_PATTYN = 4
  type(option_choice), dimension(5), parameter :: sigma_method_choices = [ &
       option_choice(SIGMA_METHOD_GLIDE,   'glide',   'standard Glide sigma levels'), &
       option_choice(SIGMA_METHOD_FILE,    'file',    'sigma levels from sigma_file'), &
       option_choice(SIGMA_METHOD_CONFIG,  'config',  'sigma levels from [sigma]'), &
       option_choice(SIGMA_METHOD_UNIFORM, 'uniform', 'evenly spaced sigma levels'), &
       option_choice(SIGMA_METHOD_PATTYN,  'pattyn',  'Pattyn sigma levels') ]

  ! [time] dt_input_option
  integer, parameter :: DT_INPUT_OPTION_DT = 0
  integer, parameter :: DT_INPUT_OPTION_NSTEPS_PER_YEAR = 1
  type(option_choice), dimension(2), parameter :: dt_input_option_choices = [ &
       option_choice(DT_INPUT_OPTION_DT,              'dt',              'time step set by dt (yr)'), &
       option_choice(DT_INPUT_OPTION_NSTEPS_PER_YEAR, 'nsteps_per_year', 'time step set by nsteps_per_year') ]

  ! [diagnostics] dm_dt_diag
  integer, parameter :: DM_DT_DIAG_KG_PER_S = 0
  integer, parameter :: DM_DT_DIAG_GT_PER_YR = 1
  type(option_choice), dimension(2), parameter :: dm_dt_diag_choices = [ &
       option_choice(DM_DT_DIAG_KG_PER_S,  'kg_per_s',  'write dmass/dt diagnostic in units of kg/s'), &
       option_choice(DM_DT_DIAG_GT_PER_YR, 'gt_per_yr', 'write dmass/dt diagnostic in units of Gt/yr') ]

  ! [diagnostics] diag_ice_threshold
  integer, parameter :: DIAG_ICE_THRESHOLD_ZERO = 0
  integer, parameter :: DIAG_ICE_THRESHOLD_THCKLIM = 1
  type(option_choice), dimension(2), parameter :: diag_ice_threshold_choices = [ &
       option_choice(DIAG_ICE_THRESHOLD_ZERO,    'zero',    'include cells with H > 0 in global diagnostics'), &
       option_choice(DIAG_ICE_THRESHOLD_THCKLIM, 'thcklim', 'include cells with H > thcklim in global diagnostics') ]

  ! [options] dycore
  ! Read by GetDycore (glimmer_config), which also accepts the integers; the choices are used for the log.
  ! Note: DYCORE_GLIDE is not a Glissade option, but it is used in ncdf_template.F90.in,
  !       from which both glide_io and glissade_io are generated.
  integer, parameter :: DYCORE_GLIDE = 0        ! old shallow-ice dycore from Glimmer
  integer, parameter :: DYCORE_GLISSADE = 2     ! finite-element solver
  type(option_choice), dimension(1), parameter :: dycore_choices = [ &
       option_choice(DYCORE_GLISSADE, 'glissade', '') ]

  ! [options] restart
  ! The constants NO_RESTART, STANDARD_RESTART and HYBRID_RESTART (glimmer_ncdf) are shared with Glide.
  type(option_choice), dimension(3), parameter :: restart_choices = [ &
       option_choice(NO_RESTART,       'none',     'not a restart'), &
       option_choice(STANDARD_RESTART, 'standard', 'standard restart'), &
       option_choice(HYBRID_RESTART,   'hybrid',   'hybrid restart (model time is reset)') ]


  ! [numerics] compute_blocks
  integer, parameter :: COMPUTE_BLOCKS_ALL = 0          ! compute on all blocks in the global domain
  integer, parameter :: COMPUTE_BLOCKS_ACTIVE_ONLY = 1  ! compute on active blocks only (where ice is potentially present)
  integer, parameter :: COMPUTE_BLOCKS_INQUIRE = 2      ! inquire how many active blocks there are for a given ice domain mask and block size
  type(option_choice), dimension(3), parameter :: compute_blocks_choices = [ &
       option_choice(COMPUTE_BLOCKS_ALL,         'all',         'compute on all blocks'), &
       option_choice(COMPUTE_BLOCKS_ACTIVE_ONLY, 'active_only', 'compute on active blocks only'), &
       option_choice(COMPUTE_BLOCKS_INQUIRE,     'inquire',     'inquire the number of active blocks') ]

  ! [velocity] stokes_approx
  integer, parameter :: STOKES_APPROX_LOCAL_SIA = 0
  integer, parameter :: STOKES_APPROX_SIA = 1
  integer, parameter :: STOKES_APPROX_SSA = 2
  integer, parameter :: STOKES_APPROX_BLATTER_PATTYN = 3
  integer, parameter :: STOKES_APPROX_L1L2 = 4
  integer, parameter :: STOKES_APPROX_DIVA = 5
  integer, parameter :: STOKES_APPROX_HYBRID = 6
  type(option_choice), dimension(7), parameter :: stokes_approx_choices = [ &
       option_choice(STOKES_APPROX_LOCAL_SIA,      'local_sia',      'SIA only (glissade_velo_sia)'), &
       option_choice(STOKES_APPROX_SIA,            'sia',            'SIA only (glissade_velo_higher)'), &
       option_choice(STOKES_APPROX_SSA,            'ssa',            'SSA only (glissade_velo_higher)'), &
       option_choice(STOKES_APPROX_BLATTER_PATTYN, 'blatter_pattyn', 'Blatter-Pattyn HO (glissade_velo_higher)'), &
       option_choice(STOKES_APPROX_L1L2,           'l1l2',           'depth-integrated L1L2 (glissade_velo_higher)'), &
       option_choice(STOKES_APPROX_DIVA,           'diva',           'depth-integrated viscosity (glissade_velo_higher)'), &
       option_choice(STOKES_APPROX_HYBRID,         'hybrid',         'hybrid SIA/SSA') ]

  ! [velocity] efvs_method
  integer, parameter :: EFVS_METHOD_CONSTANT = 0
  integer, parameter :: EFVS_METHOD_FLOW_FACTOR = 1
  integer, parameter :: EFVS_METHOD_NONLINEAR = 2
  type(option_choice), dimension(3), parameter :: efvs_method_choices = [ &
       option_choice(EFVS_METHOD_CONSTANT,    'constant',    'constant value'), &
       option_choice(EFVS_METHOD_FLOW_FACTOR, 'flow_factor', 'multiple of flow factor'), &
       option_choice(EFVS_METHOD_NONLINEAR,   'nonlinear',   'nonlinear, from effective strain rate') ]

  ! [velocity] surface_gradient_method
  integer, parameter :: SURFACE_GRADIENT_METHOD_CENTERED = 0
  integer, parameter :: SURFACE_GRADIENT_METHOD_UPSTREAM1 = 1
  integer, parameter :: SURFACE_GRADIENT_METHOD_UPSTREAM2 = 2
  type(option_choice), dimension(3), parameter :: surface_gradient_method_choices = [ &
       option_choice(SURFACE_GRADIENT_METHOD_CENTERED,  'centered',  'centered gradient'), &
       option_choice(SURFACE_GRADIENT_METHOD_UPSTREAM1, 'upstream1', 'first-order upstream gradient'), &
       option_choice(SURFACE_GRADIENT_METHOD_UPSTREAM2, 'upstream2', 'second-order upstream gradient') ]

  ! [velocity] surface_gradient_margin_method
  integer, parameter :: SURFACE_GRADIENT_MARGIN_METHOD_LAND = 0
  integer, parameter :: SURFACE_GRADIENT_MARGIN_METHOD_HYBRID = 1
  integer, parameter :: SURFACE_GRADIENT_MARGIN_METHOD_MARINE = 2
  type(option_choice), dimension(3), parameter :: surface_gradient_margin_method_choices = [ &
       option_choice(SURFACE_GRADIENT_MARGIN_METHOD_LAND,   'land',   'compute edge gradient when either cell has ice'), &
       option_choice(SURFACE_GRADIENT_MARGIN_METHOD_HYBRID, 'hybrid', 'compute edge gradient when ice lies above ice-free land'), &
       option_choice(SURFACE_GRADIENT_MARGIN_METHOD_MARINE, 'marine', 'compute edge gradient when both cells have ice') ]

  ! [velocity_solver] nonlinear_method
  integer, parameter :: NONLINEAR_METHOD_PICARD = 0
  integer, parameter :: NONLINEAR_METHOD_ACCELERATED_PICARD = 1
  type(option_choice), dimension(2), parameter :: nonlinear_method_choices = [ &
       option_choice(NONLINEAR_METHOD_PICARD,             'picard',             'standard Picard iteration'), &
       option_choice(NONLINEAR_METHOD_ACCELERATED_PICARD, 'accelerated_picard', 'Picard iteration with acceleration') ]

  ! [velocity_solver] nonlinear_residual
  integer, parameter :: NONLINEAR_RESIDUAL_MAX = 0
  integer, parameter :: NONLINEAR_RESIDUAL_MAX_NO_UBAS = 1
  integer, parameter :: NONLINEAR_RESIDUAL_MEAN = 2
  integer, parameter :: NONLINEAR_RESIDUAL_L2NORM = 3
  integer, parameter :: NONLINEAR_RESIDUAL_RELATIVE_L2NORM = 4
  type(option_choice), dimension(5), parameter :: nonlinear_residual_choices = [ &
       option_choice(NONLINEAR_RESIDUAL_MAX,             'max',             'max value'), &
       option_choice(NONLINEAR_RESIDUAL_MAX_NO_UBAS,     'max_no_ubas',     'max value ignoring ubas'), &
       option_choice(NONLINEAR_RESIDUAL_MEAN,            'mean',            'mean value'), &
       option_choice(NONLINEAR_RESIDUAL_L2NORM,          'l2norm',          'L2 norm of Ax-b = resid'), &
       option_choice(NONLINEAR_RESIDUAL_RELATIVE_L2NORM, 'relative_l2norm', 'relative L2 norm, |Ax-b|/|b|') ]

  ! [velocity_solver] linear_solver
  ! The SLAP solvers (one task only) are in glimmer_solve; glissade_velo_higher_slap maps the
  ! LINEAR_SOLVER_SLAP_* choices to that library's methods (SLAP_*).
  integer, parameter :: LINEAR_SOLVER_SLAP_PCG = 0
  integer, parameter :: LINEAR_SOLVER_SLAP_BICG = 1
  integer, parameter :: LINEAR_SOLVER_SLAP_GMRES = 2
  integer, parameter :: LINEAR_SOLVER_PCG = 3
  integer, parameter :: LINEAR_SOLVER_PCG_CHRONGEAR = 4
  type(option_choice), dimension(5), parameter :: linear_solver_choices = [ &
       option_choice(LINEAR_SOLVER_SLAP_PCG,      'slap_pcg',      'SLAP PCG with incomplete Cholesky preconditioner'), &
       option_choice(LINEAR_SOLVER_SLAP_BICG,     'slap_bicg',     'SLAP BiCG with LU preconditioner'), &
       option_choice(LINEAR_SOLVER_SLAP_GMRES,    'slap_gmres',    'SLAP GMRES with LU preconditioner'), &
       option_choice(LINEAR_SOLVER_PCG,           'pcg',           'native PCG solver, standard'), &
       option_choice(LINEAR_SOLVER_PCG_CHRONGEAR, 'pcg_chrongear', 'native PCG solver, Chronopoulos-Gear') ]

  ! [velocity_solver] preconditioner
  integer, parameter :: PRECONDITIONER_NONE = 0
  integer, parameter :: PRECONDITIONER_DIAGONAL = 1
  integer, parameter :: PRECONDITIONER_SIA  = 2
  integer, parameter :: PRECONDITIONER_TRIDIAG_LOCAL  = 3
  integer, parameter :: PRECONDITIONER_TRIDIAG_GLOBAL = 4
  type(option_choice), dimension(5), parameter :: preconditioner_choices = [ &
       option_choice(PRECONDITIONER_NONE,           'none',           'no preconditioner (native PCG)'), &
       option_choice(PRECONDITIONER_DIAGONAL,       'diagonal',       'diagonal preconditioner (native PCG)'), &
       option_choice(PRECONDITIONER_SIA,            'sia',            'SIA preconditioner (native PCG)'), &
       option_choice(PRECONDITIONER_TRIDIAG_LOCAL,  'tridiag_local',  'local tridiagonal preconditioner (native PCG)'), &
       option_choice(PRECONDITIONER_TRIDIAG_GLOBAL, 'tridiag_global', 'global tridiagonal preconditioner (native PCG)') ]

  ! [velocity_solver] assemble_beta, assemble_taud, assemble_bfric, assemble_lateral
  integer, parameter :: ASSEMBLE_BETA_STANDARD = 0
  integer, parameter :: ASSEMBLE_BETA_LOCAL = 1
  type(option_choice), dimension(2), parameter :: assemble_beta_choices = [ &
       option_choice(ASSEMBLE_BETA_STANDARD, 'standard', 'standard finite-element assembly'), &
       option_choice(ASSEMBLE_BETA_LOCAL,    'local',    'use local beta at each vertex') ]

  integer, parameter :: ASSEMBLE_TAUD_STANDARD = 0
  integer, parameter :: ASSEMBLE_TAUD_LOCAL = 1
  type(option_choice), dimension(2), parameter :: assemble_taud_choices = [ &
       option_choice(ASSEMBLE_TAUD_STANDARD, 'standard', 'standard finite-element assembly'), &
       option_choice(ASSEMBLE_TAUD_LOCAL,    'local',    'use local driving stress at each vertex') ]

  integer, parameter :: ASSEMBLE_BFRIC_STANDARD = 0
  integer, parameter :: ASSEMBLE_BFRIC_LOCAL = 1
  type(option_choice), dimension(2), parameter :: assemble_bfric_choices = [ &
       option_choice(ASSEMBLE_BFRIC_STANDARD, 'standard', 'standard finite-element assembly'), &
       option_choice(ASSEMBLE_BFRIC_LOCAL,    'local',    'use local basal friction at each vertex') ]

  integer, parameter :: ASSEMBLE_LATERAL_STANDARD = 0
  integer, parameter :: ASSEMBLE_LATERAL_LOCAL = 1
  type(option_choice), dimension(2), parameter :: assemble_lateral_choices = [ &
       option_choice(ASSEMBLE_LATERAL_STANDARD, 'standard', 'standard finite-element assembly'), &
       option_choice(ASSEMBLE_LATERAL_LOCAL,    'local',    'use local thck and usrf on each cell face') ]

  ! [advection] advection_method (values 0-2 were Glide-only)
  integer, parameter :: ADVECTION_METHOD_REMAP = 3
  integer, parameter :: ADVECTION_METHOD_UPWIND = 4
  integer, parameter :: ADVECTION_METHOD_FIXED_THICKNESS = 5
  type(option_choice), dimension(3), parameter :: advection_method_choices = [ &
       option_choice(ADVECTION_METHOD_REMAP,           'remap',           'incremental remapping'), &
       option_choice(ADVECTION_METHOD_UPWIND,          'upwind',          'first-order upwind'), &
       option_choice(ADVECTION_METHOD_FIXED_THICKNESS, 'fixed_thickness', 'thickness fixed at initial value') ]

  ! [advection] vertical_remap_order
  integer, parameter :: VERTICAL_REMAP_ORDER_FIRST = 0
  integer, parameter :: VERTICAL_REMAP_ORDER_SECOND = 1
  type(option_choice), dimension(2), parameter :: vertical_remap_order_choices = [ &
       option_choice(VERTICAL_REMAP_ORDER_FIRST,  'first',  'first-order accurate'), &
       option_choice(VERTICAL_REMAP_ORDER_SECOND, 'second', 'second-order accurate') ]


  ! [temperature] temperature_method (enthalpy is not scientifically supported)
  integer, parameter :: TEMPERATURE_METHOD_ISOTHERMAL_ARTM = 0
  integer, parameter :: TEMPERATURE_METHOD_PROGNOSTIC = 1
  integer, parameter :: TEMPERATURE_METHOD_FIXED = 2
  integer, parameter :: TEMPERATURE_METHOD_ENTHALPY = 3
  type(option_choice), dimension(4), parameter :: temperature_method_choices = [ &
       option_choice(TEMPERATURE_METHOD_ISOTHERMAL_ARTM, 'isothermal_artm', 'isothermal, set to surface air temperature'), &
       option_choice(TEMPERATURE_METHOD_PROGNOSTIC,      'prognostic',      'prognostic temperature'), &
       option_choice(TEMPERATURE_METHOD_FIXED,           'fixed',           'constant in time'), &
       option_choice(TEMPERATURE_METHOD_ENTHALPY,        'enthalpy',        'prognostic enthalpy') ]

  ! [temperature] temperature_init
  integer, parameter :: TEMPERATURE_INIT_ZERO_CELSIUS = 0
  integer, parameter :: TEMPERATURE_INIT_ARTM = 1
  integer, parameter :: TEMPERATURE_INIT_LINEAR = 2
  integer, parameter :: TEMPERATURE_INIT_ADVECTIVE_DIFFUSIVE = 3
  integer, parameter :: TEMPERATURE_INIT_FILE = 4
  type(option_choice), dimension(5), parameter :: temperature_init_choices = [ &
       option_choice(TEMPERATURE_INIT_ZERO_CELSIUS,        'zero_celsius',        'set to 0 C'), &
       option_choice(TEMPERATURE_INIT_ARTM,                'artm',                'set to surface air temperature'), &
       option_choice(TEMPERATURE_INIT_LINEAR,              'linear',              'linear vertical profile'), &
       option_choice(TEMPERATURE_INIT_ADVECTIVE_DIFFUSIVE, 'advective_diffusive', 'advective-diffusive balance'), &
       option_choice(TEMPERATURE_INIT_FILE,                'file',                'temperature from an external file') ]

  ! [temperature] dissipation_method
  integer, parameter :: DISSIPATION_METHOD_NONE = 0
  integer, parameter :: DISSIPATION_METHOD_SIA = 1
  integer, parameter :: DISSIPATION_METHOD_FIRST_ORDER = 2
  type(option_choice), dimension(3), parameter :: dissipation_method_choices = [ &
       option_choice(DISSIPATION_METHOD_NONE,        'none',        'no dissipation'), &
       option_choice(DISSIPATION_METHOD_SIA,         'sia',         '0-order SIA'), &
       option_choice(DISSIPATION_METHOD_FIRST_ORDER, 'first_order', 'first-order model (Blatter-Pattyn)') ]

  ! [flow_law] flow_law
  integer, parameter :: FLOW_LAW_CONSTANT = 0
  integer, parameter :: FLOW_LAW_PATERSON_BUDD_MINUS5C = 1
  integer, parameter :: FLOW_LAW_PATERSON_BUDD = 2
  type(option_choice), dimension(3), parameter :: flow_law_choices = [ &
       option_choice(FLOW_LAW_CONSTANT,              'constant',              'constant flow factor flwa_const'), &
       option_choice(FLOW_LAW_PATERSON_BUDD_MINUS5C, 'paterson_budd_minus5c', 'Paterson and Budd (T = -5 C)'), &
       option_choice(FLOW_LAW_PATERSON_BUDD,         'paterson_budd',         'Paterson and Budd') ]

  ! [basal_physics] sliding_law
  integer, parameter :: SLIDING_LAW_BETA_CONSTANT = 0
  integer, parameter :: SLIDING_LAW_BETA_BPMP = 1
  integer, parameter :: SLIDING_LAW_PSEUDO_PLASTIC = 2
  integer, parameter :: SLIDING_LAW_PSEUDO_PLASTIC_OLD = 3
  integer, parameter :: SLIDING_LAW_NO_SLIP = 4
  integer, parameter :: SLIDING_LAW_BETA_FILE = 5
  integer, parameter :: SLIDING_LAW_NO_SLIP_DIRICHLET = 6
  integer, parameter :: SLIDING_LAW_ZOET_IVERSON = 7
  integer, parameter :: SLIDING_LAW_ISMIP_HOM_C = 8
  integer, parameter :: SLIDING_LAW_POWERLAW = 9
  integer, parameter :: SLIDING_LAW_COULOMB = 10
  integer, parameter :: SLIDING_LAW_SCHOOF = 11
  integer, parameter :: SLIDING_LAW_TSAI = 12
  integer, parameter :: SLIDING_LAW_POWERLAW_EFFECPRESS = 13
  integer, parameter :: SLIDING_LAW_YIELD_PICARD = 14
  type(option_choice), dimension(15), parameter :: sliding_law_choices = [ &
       option_choice(SLIDING_LAW_BETA_CONSTANT,       'beta_constant',       'constant beta'), &
       option_choice(SLIDING_LAW_BETA_BPMP,           'beta_bpmp',           'beta depends on basal temperature (melting or frozen)'), &
       option_choice(SLIDING_LAW_PSEUDO_PLASTIC,      'pseudo_plastic',      'pseudo-plastic sliding law, new C_c options'), &
       option_choice(SLIDING_LAW_PSEUDO_PLASTIC_OLD,  'pseudo_plastic_old',  'pseudo-plastic sliding law, old tan(phi) options'), &
       option_choice(SLIDING_LAW_NO_SLIP,             'no_slip',             'no slip (using large beta)'), &
       option_choice(SLIDING_LAW_BETA_FILE,           'beta_file',           'beta from external file'), &
       option_choice(SLIDING_LAW_NO_SLIP_DIRICHLET,   'no_slip_dirichlet',   'no slip (Dirichlet implementation)'), &
       option_choice(SLIDING_LAW_ZOET_IVERSON,        'zoet_iverson',        'Zoet-Iverson sliding law'), &
       option_choice(SLIDING_LAW_ISMIP_HOM_C,         'ismip_hom_c',         'beta as in ISMIP-HOM test C'), &
       option_choice(SLIDING_LAW_POWERLAW,            'powerlaw',            'power law'), &
       option_choice(SLIDING_LAW_COULOMB,             'coulomb',             'Coulomb friction law with effective pressure'), &
       option_choice(SLIDING_LAW_SCHOOF,              'schoof',              'Schoof friction law'), &
       option_choice(SLIDING_LAW_TSAI,                'tsai',                'min of Coulomb stress and power-law stress (Tsai)'), &
       option_choice(SLIDING_LAW_POWERLAW_EFFECPRESS, 'powerlaw_effecpress', 'power law using effective pressure'), &
       option_choice(SLIDING_LAW_YIELD_PICARD,        'yield_picard',        'till yield stress (Picard)') ]

  ! [basal_physics] btrc_method (local SIA solver only)
  integer, parameter :: BTRC_METHOD_NONE = 0
  integer, parameter :: BTRC_METHOD_CONSTANT = 1
  integer, parameter :: BTRC_METHOD_CONSTANT_BWAT = 2
  integer, parameter :: BTRC_METHOD_CONSTANT_BPMP = 3
  type(option_choice), dimension(4), parameter :: btrc_method_choices = [ &
       option_choice(BTRC_METHOD_NONE,          'none',          'no basal sliding'), &
       option_choice(BTRC_METHOD_CONSTANT,      'constant',      'constant basal traction'), &
       option_choice(BTRC_METHOD_CONSTANT_BWAT, 'constant_bwat', 'constant where bwat > 0'), &
       option_choice(BTRC_METHOD_CONSTANT_BPMP, 'constant_bpmp', 'constant where T = Tpmp') ]

  ! [basal_physics] beta_limit_method
  integer, parameter :: BETA_LIMIT_METHOD_ABSOLUTE = 0
  integer, parameter :: BETA_LIMIT_METHOD_FGROUND_SCALED = 1
  type(option_choice), dimension(2), parameter :: beta_limit_method_choices = [ &
       option_choice(BETA_LIMIT_METHOD_ABSOLUTE,       'absolute',       'absolute beta limit based on beta_grounded_min'), &
       option_choice(BETA_LIMIT_METHOD_FGROUND_SCALED, 'fground_scaled', 'beta is limited, then scaled by f_ground_cell') ]

  ! [basal_physics] powerlaw_c_method
  integer, parameter :: POWERLAW_C_METHOD_CONSTANT = 0
  integer, parameter :: POWERLAW_C_METHOD_INVERSION = 1
  integer, parameter :: POWERLAW_C_METHOD_FILE = 2
  integer, parameter :: POWERLAW_C_METHOD_INVERSION_BASIN = 3
  integer, parameter :: POWERLAW_C_METHOD_FROM_COULOMB_C = 4
  type(option_choice), dimension(5), parameter :: powerlaw_c_method_choices = [ &
       option_choice(POWERLAW_C_METHOD_CONSTANT,        'constant',        'spatially uniform friction parameter Cp'), &
       option_choice(POWERLAW_C_METHOD_INVERSION,       'inversion',       'invert for 2D friction parameter Cp'), &
       option_choice(POWERLAW_C_METHOD_FILE,            'file',            'read friction parameter Cp from file'), &
       option_choice(POWERLAW_C_METHOD_INVERSION_BASIN, 'inversion_basin', 'invert for basin-scale friction parameter Cp'), &
       option_choice(POWERLAW_C_METHOD_FROM_COULOMB_C,  'from_coulomb_c',  'friction parameter Cp is a function of Cc') ]

  ! [basal_physics] coulomb_c_method
  integer, parameter :: COULOMB_C_METHOD_CONSTANT = 0
  integer, parameter :: COULOMB_C_METHOD_INVERSION = 1
  integer, parameter :: COULOMB_C_METHOD_FILE = 2
  integer, parameter :: COULOMB_C_METHOD_INVERSION_BASIN = 3
  type(option_choice), dimension(4), parameter :: coulomb_c_method_choices = [ &
       option_choice(COULOMB_C_METHOD_CONSTANT,        'constant',        'spatially uniform friction parameter Cc'), &
       option_choice(COULOMB_C_METHOD_INVERSION,       'inversion',       'invert for 2D friction parameter Cc'), &
       option_choice(COULOMB_C_METHOD_FILE,            'file',            'read friction parameter Cc from file'), &
       option_choice(COULOMB_C_METHOD_INVERSION_BASIN, 'inversion_basin', 'invert for basin-scale friction parameter Cc') ]

  ! [basal_physics] effecpress_method
  integer, parameter :: EFFECPRESS_METHOD_OVERBURDEN = 0
  integer, parameter :: EFFECPRESS_METHOD_BPMP = 1
  integer, parameter :: EFFECPRESS_METHOD_BWAT = 2
  integer, parameter :: EFFECPRESS_METHOD_CAVITY_SHEET = 3
  integer, parameter :: EFFECPRESS_METHOD_BWAT_BVP = 4
  type(option_choice), dimension(5), parameter :: effecpress_method_choices = [ &
       option_choice(EFFECPRESS_METHOD_OVERBURDEN,   'overburden',   'full overburden pressure'), &
       option_choice(EFFECPRESS_METHOD_BPMP,         'bpmp',         'reduced effecpress near pressure melting point'), &
       option_choice(EFFECPRESS_METHOD_BWAT,         'bwat',         'macroporous sheet; effecpress depends on bwat'), &
       option_choice(EFFECPRESS_METHOD_CAVITY_SHEET, 'cavity_sheet', 'cavity-sheet; effecpress depends on opening, closing'), &
       option_choice(EFFECPRESS_METHOD_BWAT_BVP,     'bwat_bvp',     'reduced effecpress where bwat > 0 (B/vP)') ]

  ! [basal_hydro] bwat_method
  integer, parameter :: BWAT_METHOD_NONE = 0
  integer, parameter :: BWAT_METHOD_CONSTANT = 1
  integer, parameter :: BWAT_METHOD_LOCAL_TILL = 2
  integer, parameter :: BWAT_METHOD_FLUX_ROUTING = 3
  type(option_choice), dimension(4), parameter :: bwat_method_choices = [ &
       option_choice(BWAT_METHOD_NONE,         'none',         'zero basal water depth'), &
       option_choice(BWAT_METHOD_CONSTANT,     'constant',     'constant basal water depth'), &
       option_choice(BWAT_METHOD_LOCAL_TILL,   'local_till',   'basal water depth computed from local till model'), &
       option_choice(BWAT_METHOD_FLUX_ROUTING, 'flux_routing', 'steady-state flux routing') ]

  ! [basal_hydro] flux_routing_method
  integer, parameter :: FLUX_ROUTING_METHOD_D8 = 0
  integer, parameter :: FLUX_ROUTING_METHOD_DINF = 1
  integer, parameter :: FLUX_ROUTING_METHOD_FD8 = 2
  type(option_choice), dimension(3), parameter :: flux_routing_method_choices = [ &
       option_choice(FLUX_ROUTING_METHOD_D8,   'd8',   'D8; route flux to lowest-elevation neighbor'), &
       option_choice(FLUX_ROUTING_METHOD_DINF, 'dinf', 'Dinf; route flux to two lower-elevation neighbors'), &
       option_choice(FLUX_ROUTING_METHOD_FD8,  'fd8',  'FD8; route flux to all lower-elevation neighbors') ]

  ! [basal_hydro] cavity_open_slide (only with flux routing and the cavity-sheet model)
  integer, parameter :: CAVITY_OPEN_SLIDE_NONE = 0
  integer, parameter :: CAVITY_OPEN_SLIDE_FIXED_UB = 1
  integer, parameter :: CAVITY_OPEN_SLIDE_DYNAMIC_UB = 2
  type(option_choice), dimension(3), parameter :: cavity_open_slide_choices = [ &
       option_choice(CAVITY_OPEN_SLIDE_NONE,       'none',       'no opening by sliding'), &
       option_choice(CAVITY_OPEN_SLIDE_FIXED_UB,   'fixed_ub',   'opening by sliding with fixed u_b'), &
       option_choice(CAVITY_OPEN_SLIDE_DYNAMIC_UB, 'dynamic_ub', 'opening by sliding with dynamic u_b') ]

  ! [basal_hydro] cavity_open_melt (only with flux routing and the cavity-sheet model)
  integer, parameter :: CAVITY_OPEN_MELT_NONE = 0
  integer, parameter :: CAVITY_OPEN_MELT_BMLT = 1
  integer, parameter :: CAVITY_OPEN_MELT_DISSIPATION = 2
  integer, parameter :: CAVITY_OPEN_MELT_BMLT_DISSIPATION = 3
  type(option_choice), dimension(4), parameter :: cavity_open_melt_choices = [ &
       option_choice(CAVITY_OPEN_MELT_NONE,             'none',             'no opening by melting'), &
       option_choice(CAVITY_OPEN_MELT_BMLT,             'bmlt',             'opening by melting based on bmlt_ground'), &
       option_choice(CAVITY_OPEN_MELT_DISSIPATION,      'dissipation',      'opening by melting based on cavity dissipation'), &
       option_choice(CAVITY_OPEN_MELT_BMLT_DISSIPATION, 'bmlt_dissipation', 'opening based on bmlt_ground plus dissipation') ]

  ! [grounding_line] grounding_line_param
  integer, parameter :: GROUNDING_LINE_PARAM_NONE = 0
  integer, parameter :: GROUNDING_LINE_PARAM_VERTICES = 1
  integer, parameter :: GROUNDING_LINE_PARAM_DELUXE = 2
  type(option_choice), dimension(3), parameter :: grounding_line_param_choices = [ &
       option_choice(GROUNDING_LINE_PARAM_NONE,     'none',     'f_ground = 0 or 1; no GLP'), &
       option_choice(GROUNDING_LINE_PARAM_VERTICES, 'vertices', 'GLP for basal friction with 0 <= f_ground <= 1 for vertices'), &
       option_choice(GROUNDING_LINE_PARAM_DELUXE,   'deluxe',   'deluxe GLP, 0 <= f_ground <= 1 for both vertices and cells') ]

  ! [grounding_line] fground_no_glp
  integer, parameter :: FGROUND_NO_GLP_GROUNDED_NEIGHBOR = 0
  integer, parameter :: FGROUND_NO_GLP_FLOTATION_FUNCTION = 1
  type(option_choice), dimension(2), parameter :: fground_no_glp_choices = [ &
       option_choice(FGROUND_NO_GLP_GROUNDED_NEIGHBOR,  'grounded_neighbor',  'f_ground = 1 at vertex if any neighbor cell is grounded'), &
       option_choice(FGROUND_NO_GLP_FLOTATION_FUNCTION, 'flotation_function', 'f_ground = 0 or 1 at vertex based on staggered f_flotation') ]

  ! [grounding_line] grounding_line_bmlt
  integer, parameter :: GROUNDING_LINE_BMLT_NONE = 0
  integer, parameter :: GROUNDING_LINE_BMLT_FLOATING_FRACTION = 1
  integer, parameter :: GROUNDING_LINE_BMLT_ZERO_GROUNDED = 2
  type(option_choice), dimension(3), parameter :: grounding_line_bmlt_choices = [ &
       option_choice(GROUNDING_LINE_BMLT_NONE,              'none',              'no GLP for bmlt_float'), &
       option_choice(GROUNDING_LINE_BMLT_FLOATING_FRACTION, 'floating_fraction', 'weigh bmlt_float by floating fraction of cell'), &
       option_choice(GROUNDING_LINE_BMLT_ZERO_GROUNDED,     'zero_grounded',     'set bmlt_float = 0 in partly grounded cells') ]

  ! [grounding_line] flotation_function
  integer, parameter :: FLOTATION_FUNCTION_PATTYN = 0
  integer, parameter :: FLOTATION_FUNCTION_INVERSE_PATTYN = 1
  integer, parameter :: FLOTATION_FUNCTION_LINEAR = 2
  type(option_choice), dimension(3), parameter :: flotation_function_choices = [ &
       option_choice(FLOTATION_FUNCTION_PATTYN,         'pattyn',         'f_pattyn = (-rhoo*b)/(rhoi*H)'), &
       option_choice(FLOTATION_FUNCTION_INVERSE_PATTYN, 'inverse_pattyn', '1/fpattyn = (rhoi*H)/(-rhoo*b)'), &
       option_choice(FLOTATION_FUNCTION_LINEAR,         'linear',         'linear = -b - (rhoi/rhoo)*H') ]

  ! [basal_melt] bmlt_float_method (value 5 was the old MISOMIP plume model)
  integer, parameter :: BMLT_FLOAT_METHOD_NONE = 0
  integer, parameter :: BMLT_FLOAT_METHOD_MISMIP = 1
  integer, parameter :: BMLT_FLOAT_METHOD_CONSTANT = 2
  integer, parameter :: BMLT_FLOAT_METHOD_DEPTH = 3
  integer, parameter :: BMLT_FLOAT_METHOD_FILE = 4
  integer, parameter :: BMLT_FLOAT_METHOD_THERMAL_FORCING = 5
  type(option_choice), dimension(6), parameter :: bmlt_float_method_choices = [ &
       option_choice(BMLT_FLOAT_METHOD_NONE,            'none',            'no melting of floating ice'), &
       option_choice(BMLT_FLOAT_METHOD_MISMIP,          'mismip',          'MISMIP+ melt rate profile'), &
       option_choice(BMLT_FLOAT_METHOD_CONSTANT,        'constant',        'constant melt rate'), &
       option_choice(BMLT_FLOAT_METHOD_DEPTH,           'depth',           'depth-dependent melt rate'), &
       option_choice(BMLT_FLOAT_METHOD_FILE,            'file',            'melt rate from external file'), &
       option_choice(BMLT_FLOAT_METHOD_THERMAL_FORCING, 'thermal_forcing', 'melt rate from thermal forcing') ]

  ! [basal_melt] bmlt_float_thermal_forcing_param
  integer, parameter :: BMLT_FLOAT_TF_PARAM_QUADRATIC = 0
  integer, parameter :: BMLT_FLOAT_TF_PARAM_ISMIP6_LOCAL = 1
  integer, parameter :: BMLT_FLOAT_TF_PARAM_ISMIP6_NONLOCAL = 2
  integer, parameter :: BMLT_FLOAT_TF_PARAM_ISMIP6_NONLOCAL_SLOPE = 3
  type(option_choice), dimension(4), parameter :: bmlt_float_thermal_forcing_param_choices = [ &
       option_choice(BMLT_FLOAT_TF_PARAM_QUADRATIC,             'quadratic',             'quadratic function of thermal forcing'), &
       option_choice(BMLT_FLOAT_TF_PARAM_ISMIP6_LOCAL,          'ismip6_local',          'ISMIP6 local quadratic'), &
       option_choice(BMLT_FLOAT_TF_PARAM_ISMIP6_NONLOCAL,       'ismip6_nonlocal',       'ISMIP6 nonlocal quadratic'), &
       option_choice(BMLT_FLOAT_TF_PARAM_ISMIP6_NONLOCAL_SLOPE, 'ismip6_nonlocal_slope', 'ISMIP6 nonlocal quadratic, slope-dependent') ]

  ! [basal_melt] deltaT_ocn_method
  integer, parameter :: DELTAT_OCN_METHOD_NONE = 0
  integer, parameter :: DELTAT_OCN_METHOD_INVERSION = 1
  integer, parameter :: DELTAT_OCN_METHOD_FILE = 2
  integer, parameter :: DELTAT_OCN_METHOD_INVERSION_BASIN = 3
  integer, parameter :: DELTAT_OCN_METHOD_CALIBRATE_BASIN = 4
  integer, parameter :: DELTAT_OCN_METHOD_DTHCK_DT = 5
  type(option_choice), dimension(6), parameter :: deltaT_ocn_method_choices = [ &
       option_choice(DELTAT_OCN_METHOD_NONE,            'none',            'deltaT_ocn = 0'), &
       option_choice(DELTAT_OCN_METHOD_INVERSION,       'inversion',       'invert for 2D deltaT_ocn based on thck'), &
       option_choice(DELTAT_OCN_METHOD_FILE,            'file',            'read deltaT_ocn from external file'), &
       option_choice(DELTAT_OCN_METHOD_INVERSION_BASIN, 'inversion_basin', 'invert for basin-scale deltaT_ocn'), &
       option_choice(DELTAT_OCN_METHOD_CALIBRATE_BASIN, 'calibrate_basin', 'calibrate deltaT_basin to match melt target'), &
       option_choice(DELTAT_OCN_METHOD_DTHCK_DT,        'dthck_dt',        'set deltaT_ocn based on observed dthck_dt') ]

  ! [ocean_data] ocean_data_domain
  integer, parameter :: OCEAN_DATA_DOMAIN_INTERNAL = 0
  integer, parameter :: OCEAN_DATA_DOMAIN_FILE = 1
  integer, parameter :: OCEAN_DATA_DOMAIN_COUPLER = 2
  type(option_choice), dimension(3), parameter :: ocean_data_domain_choices = [ &
       option_choice(OCEAN_DATA_DOMAIN_INTERNAL, 'internal', 'ocean data computed internally by CISM'), &
       option_choice(OCEAN_DATA_DOMAIN_FILE,     'file',     'ocean data read from external file'), &
       option_choice(OCEAN_DATA_DOMAIN_COUPLER,  'coupler',  'ocean data from coupler via Glad') ]

  ! [calving] calving_law (the old value 2 (a fixed fraction of floating ice) was removed)
  integer, parameter :: CALVING_LAW_NONE = 0
  integer, parameter :: CALVING_LAW_FLOAT_ZERO = 1
  integer, parameter :: CALVING_LAW_RELX_THRESHOLD = 2
  integer, parameter :: CALVING_LAW_TOPG_THRESHOLD = 3
  integer, parameter :: CALVING_LAW_GRID_MASK = 4
  integer, parameter :: CALVING_LAW_ADVANCE_RETREAT_RATE = 5
  integer, parameter :: CALVING_LAW_THCK_THRESHOLD = 6
  integer, parameter :: CALVING_LAW_STRESS = 7
  integer, parameter :: CALVING_LAW_STRESS_STOCHASTIC = 8
  integer, parameter :: CALVING_LAW_DAMAGE = 9
  integer, parameter :: CALVING_LAW_EIGENCALVING = 10
  type(option_choice), dimension(11), parameter :: calving_law_choices = [ &
       option_choice(CALVING_LAW_NONE,                 'none',                 'no calving law'), &
       option_choice(CALVING_LAW_FLOAT_ZERO,           'float_zero',           'remove all floating ice'), &
       option_choice(CALVING_LAW_RELX_THRESHOLD,       'relx_threshold',       'relaxed bedrock threshold'), &
       option_choice(CALVING_LAW_TOPG_THRESHOLD,       'topg_threshold',       'present bedrock threshold'), &
       option_choice(CALVING_LAW_GRID_MASK,            'grid_mask',            'calving based on grid location'), &
       option_choice(CALVING_LAW_ADVANCE_RETREAT_RATE, 'advance_retreat_rate', 'prescribe CF advance/retreat rate'), &
       option_choice(CALVING_LAW_THCK_THRESHOLD,       'thck_threshold',       'ice thickness threshold'), &
       option_choice(CALVING_LAW_STRESS,               'stress',               'deterministic stress-based calving'), &
       option_choice(CALVING_LAW_STRESS_STOCHASTIC,    'stress_stochastic',    'stochastic stress-based calving'), &
       option_choice(CALVING_LAW_DAMAGE,               'damage',               'damage-based calving scheme'), &
       option_choice(CALVING_LAW_EIGENCALVING,         'eigencalving',         'strain-rate (eigencalving)') ]

  ! [calving] calving_domain
  integer, parameter :: CALVING_DOMAIN_OCEAN_EDGE = 0
  integer, parameter :: CALVING_DOMAIN_EVERYWHERE = 1
  type(option_choice), dimension(2), parameter :: calving_domain_choices = [ &
       option_choice(CALVING_DOMAIN_OCEAN_EDGE, 'ocean_edge', 'calving only at the ocean edge'), &
       option_choice(CALVING_DOMAIN_EVERYWHERE, 'everywhere', 'calving in all cells where criterion is met') ]

  ! [calving] calvingmip_domain
  integer, parameter :: CALVINGMIP_DOMAIN_NONE = 0
  integer, parameter :: CALVINGMIP_DOMAIN_CIRCULAR = 1
  integer, parameter :: CALVINGMIP_DOMAIN_THULE = 2
  type(option_choice), dimension(3), parameter :: calvingmip_domain_choices = [ &
       option_choice(CALVINGMIP_DOMAIN_NONE,     'none',     'not a calvingMIP domain'), &
       option_choice(CALVINGMIP_DOMAIN_CIRCULAR, 'circular', 'calvingMIP circular domain'), &
       option_choice(CALVINGMIP_DOMAIN_THULE,    'thule',    'calvingMIP Thule domain') ]

  ! [lateral_melt] lateral_melt_method
  integer, parameter :: LATERAL_MELT_METHOD_NONE = 0
  integer, parameter :: LATERAL_MELT_METHOD_CONSTANT = 1
  integer, parameter :: LATERAL_MELT_METHOD_ISMIP = 2
  type(option_choice), dimension(3), parameter :: lateral_melt_method_choices = [ &
       option_choice(LATERAL_MELT_METHOD_NONE,     'none',     'no lateral melt at marine margin'), &
       option_choice(LATERAL_MELT_METHOD_CONSTANT, 'constant', 'constant lateral melt rate'), &
       option_choice(LATERAL_MELT_METHOD_ISMIP,    'ismip',    'ISMIP lateral melt from TF and discharge') ]

  ! [removal] force_retreat
  integer, parameter :: FORCE_RETREAT_NONE = 0
  integer, parameter :: FORCE_RETREAT_ALL_ICE = 1
  integer, parameter :: FORCE_RETREAT_FLOATING_ICE = 2
  type(option_choice), dimension(3), parameter :: force_retreat_choices = [ &
       option_choice(FORCE_RETREAT_NONE,         'none',         'no forced retreat'), &
       option_choice(FORCE_RETREAT_ALL_ICE,      'all_ice',      'force retreat of all ice with retreat mask'), &
       option_choice(FORCE_RETREAT_FLOATING_ICE, 'floating_ice', 'force retreat of floating ice with retreat mask') ]

  ! [climate] smb_input_units
  integer, parameter :: SMB_INPUT_UNITS_M_ICE_PER_YR = 0
  integer, parameter :: SMB_INPUT_UNITS_MM_WE_PER_YR = 1
  type(option_choice), dimension(2), parameter :: smb_input_units_choices = [ &
       option_choice(SMB_INPUT_UNITS_M_ICE_PER_YR, 'm_ice_per_yr', 'SMB input in units of m/yr ice (acab)'), &
       option_choice(SMB_INPUT_UNITS_MM_WE_PER_YR, 'mm_we_per_yr', 'SMB input in units of mm/yr w.e. (smb)') ]

  ! [climate] smb_input_function
  integer, parameter :: SMB_INPUT_FUNCTION_XY = 0
  integer, parameter :: SMB_INPUT_FUNCTION_XY_GRADZ = 1
  integer, parameter :: SMB_INPUT_FUNCTION_XYZ = 2
  integer, parameter :: SMB_INPUT_FUNCTION_PDD = 3
  type(option_choice), dimension(4), parameter :: smb_input_function_choices = [ &
       option_choice(SMB_INPUT_FUNCTION_XY,       'xy',       'SMB input as function of (x,y)'), &
       option_choice(SMB_INPUT_FUNCTION_XY_GRADZ, 'xy_gradz', 'SMB and d(SMB)/dz input as function of (x,y)'), &
       option_choice(SMB_INPUT_FUNCTION_XYZ,      'xyz',      'SMB input as function of (x,y,z)'), &
       option_choice(SMB_INPUT_FUNCTION_PDD,      'pdd',      'SMB computed with positive-degree scheme') ]

  ! [climate] artm_input_function
  integer, parameter :: ARTM_INPUT_FUNCTION_XY = 0
  integer, parameter :: ARTM_INPUT_FUNCTION_XY_GRADZ = 1
  integer, parameter :: ARTM_INPUT_FUNCTION_XYZ = 2
  integer, parameter :: ARTM_INPUT_FUNCTION_XY_LAPSE = 3
  type(option_choice), dimension(4), parameter :: artm_input_function_choices = [ &
       option_choice(ARTM_INPUT_FUNCTION_XY,       'xy',       'artm input as function of (x,y)'), &
       option_choice(ARTM_INPUT_FUNCTION_XY_GRADZ, 'xy_gradz', 'artm and d(artm)/dz input as function of (x,y)'), &
       option_choice(ARTM_INPUT_FUNCTION_XYZ,      'xyz',      'artm input as function of (x,y,z)'), &
       option_choice(ARTM_INPUT_FUNCTION_XY_LAPSE, 'xy_lapse', 'artm input as function of (x,y) with lapse rate') ]

  ! [climate] overwrite_acab
  ! (OVERWRITE_ACAB_INPUT_MASK, not OVERWRITE_ACAB_MASK: Fortran names are case-insensitive, and
  !  overwrite_acab_mask is a field)
  integer, parameter :: OVERWRITE_ACAB_NONE = 0
  integer, parameter :: OVERWRITE_ACAB_ACAB_ZERO = 1
  integer, parameter :: OVERWRITE_ACAB_THCK_MIN = 2
  integer, parameter :: OVERWRITE_ACAB_INPUT_MASK = 3
  type(option_choice), dimension(4), parameter :: overwrite_acab_choices = [ &
       option_choice(OVERWRITE_ACAB_NONE,      'none',      'do not overwrite acab anywhere'), &
       option_choice(OVERWRITE_ACAB_ACAB_ZERO, 'acab_zero', 'overwrite acab where input acab = 0'), &
       option_choice(OVERWRITE_ACAB_THCK_MIN,  'thck_min',  'overwrite acab where input thck <= overwrite_acab_thck_min'), &
       option_choice(OVERWRITE_ACAB_INPUT_MASK,      'mask',      'overwrite acab based on input mask') ]

  ! [geothermal] geothermal_method
  integer, parameter :: GEOTHERMAL_METHOD_UNIFORM = 0
  integer, parameter :: GEOTHERMAL_METHOD_FILE = 1
  type(option_choice), dimension(2), parameter :: geothermal_method_choices = [ &
       option_choice(GEOTHERMAL_METHOD_UNIFORM, 'uniform', 'uniform geothermal flux'), &
       option_choice(GEOTHERMAL_METHOD_FILE,    'file',    'read flux from file, if present') ]

  ! [isostasy] lithosphere_method
  integer, parameter :: LITHOSPHERE_METHOD_LOCAL = 0
  integer, parameter :: LITHOSPHERE_METHOD_ELASTIC = 1
  type(option_choice), dimension(2), parameter :: lithosphere_method_choices = [ &
       option_choice(LITHOSPHERE_METHOD_LOCAL,   'local',   'local lithosphere'), &
       option_choice(LITHOSPHERE_METHOD_ELASTIC, 'elastic', 'elastic lithosphere') ]

  ! [isostasy] asthenosphere_method
  integer, parameter :: ASTHENOSPHERE_METHOD_FLUID = 0
  integer, parameter :: ASTHENOSPHERE_METHOD_RELAXING = 1
  integer, parameter :: ASTHENOSPHERE_METHOD_RELAXING_LATVAR = 2
  type(option_choice), dimension(3), parameter :: asthenosphere_method_choices = [ &
       option_choice(ASTHENOSPHERE_METHOD_FLUID,           'fluid',           'fluid asthenosphere'), &
       option_choice(ASTHENOSPHERE_METHOD_RELAXING,        'relaxing',        'relaxing with constant timescale'), &
       option_choice(ASTHENOSPHERE_METHOD_RELAXING_LATVAR, 'relaxing_latvar', 'relaxing with latitude-varying timescale') ]

  ! [isostasy] input_topg_method
  integer, parameter :: INPUT_TOPG_METHOD_SEPARATE_RELX = 0
  integer, parameter :: INPUT_TOPG_METHOD_RELX_FROM_TOPG = 1
  integer, parameter :: INPUT_TOPG_METHOD_RELX_FROM_LOAD = 2
  type(option_choice), dimension(3), parameter :: input_topg_method_choices = [ &
       option_choice(INPUT_TOPG_METHOD_SEPARATE_RELX,  'separate_relx',  'read both topg and relx from input file'), &
       option_choice(INPUT_TOPG_METHOD_RELX_FROM_TOPG, 'relx_from_topg', 'set relx to input topg'), &
       option_choice(INPUT_TOPG_METHOD_RELX_FROM_LOAD, 'relx_from_load', 'compute relx = input topg + load') ]

  ! [glaciers] mu_star_method
  integer, parameter :: MU_STAR_METHOD_CONSTANT = 0
  integer, parameter :: MU_STAR_METHOD_INVERSION = 1
  integer, parameter :: MU_STAR_METHOD_FILE = 2
  type(option_choice), dimension(3), parameter :: mu_star_method_choices = [ &
       option_choice(MU_STAR_METHOD_CONSTANT,  'constant',  'spatially uniform glacier parameter mu_star'), &
       option_choice(MU_STAR_METHOD_INVERSION, 'inversion', 'glacier-specific mu_star found by inversion'), &
       option_choice(MU_STAR_METHOD_FILE,      'file',      'glacier-specific mu_star read from file') ]

  ! [glaciers] alpha_snow_method
  integer, parameter :: ALPHA_SNOW_METHOD_CONSTANT = 0
  integer, parameter :: ALPHA_SNOW_METHOD_INVERSION = 1
  integer, parameter :: ALPHA_SNOW_METHOD_FILE = 2
  type(option_choice), dimension(3), parameter :: alpha_snow_method_choices = [ &
       option_choice(ALPHA_SNOW_METHOD_CONSTANT,  'constant',  'spatially uniform glacier parameter alpha_snow'), &
       option_choice(ALPHA_SNOW_METHOD_INVERSION, 'inversion', 'glacier-specific alpha_snow found by inversion'), &
       option_choice(ALPHA_SNOW_METHOD_FILE,      'file',      'glacier-specific alpha_snow read from file') ]

  ! [glaciers] glacier_powerlaw_c_method
  integer, parameter :: GLACIER_POWERLAW_C_METHOD_CONSTANT = 0
  integer, parameter :: GLACIER_POWERLAW_C_METHOD_INVERSION = 1
  integer, parameter :: GLACIER_POWERLAW_C_METHOD_FILE = 2
  type(option_choice), dimension(3), parameter :: glacier_powerlaw_c_method_choices = [ &
       option_choice(GLACIER_POWERLAW_C_METHOD_CONSTANT,  'constant',  'spatially uniform glacier parameter Cp'), &
       option_choice(GLACIER_POWERLAW_C_METHOD_INVERSION, 'inversion', 'Cp for glaciers found by inversion'), &
       option_choice(GLACIER_POWERLAW_C_METHOD_FILE,      'file',      'Cp for glaciers read from file') ]

  ! [glaciers] snow_method
  integer, parameter :: SNOW_METHOD_SNOW = 0
  integer, parameter :: SNOW_METHOD_PRECIP_ARTM = 1
  type(option_choice), dimension(2), parameter :: snow_method_choices = [ &
       option_choice(SNOW_METHOD_SNOW,        'snow',        'read in snowfall rate directly'), &
       option_choice(SNOW_METHOD_PRECIP_ARTM, 'precip_artm', 'compute snowfall rate from precip and artm') ]


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_grid_type

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

    integer :: global_bc = GLOBAL_BC_PERIODIC

    ! added to handle the active-blocks option
    integer, dimension(:,:), pointer :: ice_domain_mask => null()  ! = 1 for cells that are potentially active

    ! mask to identify cells at the edge of the global domain
    integer, dimension(:,:), pointer :: global_edge_mask => null() ! = 1 for cells at edge of global domain

    integer :: nx_block = 0      ! user-specified block sizes
    integer :: ny_block = 0      ! one task per block; optionally, tasks not assigned to inactive blocks

    real(dp) :: dew = 20.d3      !> grid cell size in east-west direction (m)
    real(dp) :: dns = 20.d3      !> grid cell size in north-south direction (m)

    !TODO - Change default to true? Would then specify as false for idealized runs
    logical :: read_lat_lon = .false.
    !> if true, then read lat and lon fields from the input file and write to restarts

    ! Vertical coordinate ---------------------------------------------------

    integer :: sigma_method = SIGMA_METHOD_GLIDE
    !> \begin{description}
    !> \item[glide] compute standard Glimmer sigma coordinates
    !> \item[file] sigma coordinates are given in external file
    !> \item[config] sigma coordinates are given in configuration file
    !> \item[uniform] evenly spaced levels
    !> \item[pattyn] compute Pattyn sigma coordinates
    !> \end{description}

    character(len=fname_length) :: sigma_file = ''           !> file with sigma levels (sigma_method = file)

    real(dp),dimension(:),pointer :: sigma => null() !> Sigma values for vertical spacing of
                                                     !> model levels
    real(dp),dimension(:),pointer :: stagsigma => null()     !> Staggered values of sigma (layer midpts)
    real(dp),dimension(:),pointer :: stagwbndsigma => null() !> Staggered values of sigma (layer midpts) with boundaries

  end type glissade_grid_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_options_type

    !> Holds user options controlling the methods used in the ice-model integration.

    !-----------------------------------------------------------------------
    ! standard options
    !-----------------------------------------------------------------------

    integer :: dycore = DYCORE_GLISSADE

    !> Dycore (must be glissade):
    !> \begin{description} 
    !> \item[glissade] Glissade dycore (1st-order, FEM, serial (SLAP) or parallel (F90 native PCG solver) )
    !> \end{description}

    !WHL TODO - During the config reorg, Change name to smb_input_units
    logical :: forcewrite_final = .false.
    !> if true, then force a write to output and restart files when the model finishes

    integer :: restart = NO_RESTART
    !> if the run is a restart of a previous run
    !> \begin{description}
    !> \item[none] normal start-up (take init fields from .nc input file OR if absent, use default options)
    !> \item[standard] restart model from previous run (do not calc. temp, rate factor, or vel)
    !> \item[hybrid] hybrid restart; use restart from previous run as the input file, and reset the time
    !> \end{description}

    logical :: restart_extend_velo = .false.
    !> if true, write uvel_extend and vvel_extend to the restart file on the extended staggered mesh
    !> (required if restart velocities are nonzero on global boundaries);
    !> if false, write uvel and vvel on the standard staggered mesh

    !-----------------------------------------------------------------------
    ! Higher-order options associated with the glissade dycore
    !-----------------------------------------------------------------------


    !WHL TODO - Option 3 to be removed during the config reorg
    !> Flag that describes basal boundary condition for HO dyn core: 
    !> \begin{description}
    !> \item[beta_constant] spatially uniform value of beta_const
    !> \item[beta_bpmp] large value for frozen bed, lower value for bed at pressure melting point
    !> \item[pseudo_plastic] pseudo-plastic basal sliding law; new version with coulomb_c options
    !> \item[pseudo_plastic_old] pseudo-plastic basal sliding law; old version with tan(phi)
    !> \item[no_slip] very large value for beta to enforce no slip everywhere 
    !> \item[beta_file] beta field passed in from .nc input file as part of standard i/o
    !> \item[no_slip_dirichlet] no slip everywhere (using Dirichlet BC rather than large beta)
    !> \item[zoet_iverson] Zoet-Iverson law combining Coulomb and powerlaw behavior
    !> \item[ismip_hom_c] beta field as prescribed for ISMIP-HOM test C (serial only)
    !> \item[powerlaw] power law
    !> \item[coulomb] Coulomb friction law using effective pressure, with a constant basal flow factor (flwa_basal_const)
    !> \item[schoof] Schoof law that blends powerlaw and Coulomb behavior
    !> \item[tsai] basal stress is the minimum of Coulomb and power-law values, as in Tsai et al. (2015)
    !> \item[powerlaw_effecpress] power law using effective pressure
    !> \item[yield_picard] treat beta value as a till yield stress (in Pa) using Picard iteration
    !> \end{description}

    !WHL TODO - Remove this option during the config reorg?
    !           (Check whether it's needed for tests such as MISMIP)
    !TODO: Change the default to 2nd order vertical remapping
    ! WHL: Keeping 1st order vertical remapping for now, pending more testing
    !WHL TODO - Put the next few variables in a solver derived type during the config reorg
  end type glissade_options_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_geometry_type

    !> Holds fields and other information relating to the geometry of the ice sheet and bedrock.

    real(dp),dimension(:,:),pointer :: thck => null()
    !> ice thickness (m)

    real(dp),dimension(:,:),pointer :: stagthck => null()
    !> ice thickness averaged onto the staggered grid (m)

    real(dp),dimension(:,:),pointer :: usrf => null()
    ! elevation of the upper ice surface (m)

    real(dp),dimension(:,:),pointer :: lsrf => null() 
    ! elevation of the lower ice surface (m)

    real(dp),dimension(:,:),pointer :: topg => null() 
    !> elevation of the bed topography (m)

    real(dp),dimension(:,:),pointer :: usrf_obs => null()
    !> observed upper surface elevation (m)

    real(dp),dimension(:,:),pointer :: f_flotation => null() 
    !> flotation function, > 0 for floating ice and <= 0 for grounded ice

    real(dp),dimension(:,:),pointer :: f_ground => null() 
    !> The fractional area at each vertex which is grounded

    real(dp),dimension(:,:),pointer :: f_ground_cell => null() 
    !> The fractional area in each cell which is grounded

    real(dp),dimension(:,:,:),pointer :: ice_age => null()
    !> The age of a given ice layer; units are s internally, yr for I/O

    real(dp),dimension(:,:),pointer :: thck_old => null()           !> old ice thickness (m)
    real(dp),dimension(:,:),pointer :: dthck_dt => null()           !> ice thickness tendency (m/s)
    real(dp),dimension(:,:),pointer :: dthck_dt_obs => null()       !> observed rate of change of ice thickness (m/s)
    real(dp),dimension(:,:),pointer :: dthck_dt_obs_basin => null() !> basin_average of dthck_dt_obs (m/s)

    real(dp),dimension(:,:),pointer :: cell_area => null()          !> grid cell area (m^2)

    integer :: ntracers
    !> number of tracers to be transported

    real(dp),dimension(:,:,:,:), pointer :: tracers => null()
    !> all 3D tracers packed into one array

    real(dp),dimension(:,:,:), pointer :: tracers_usrf => null()
    !> value at upper surface for all tracers

    real(dp),dimension(:,:,:), pointer :: tracers_lsrf => null()
    !> value at lower surface for all tracers

    !WHL TODO - Move masks to a mask derived type (probably after the config reorg)
    integer, dimension(:,:),pointer :: ice_mask => null()        !> = 1 where ice is present, else = 0
    integer, dimension(:,:),pointer :: ice_mask_stag => null()   !> = 1 where ice is present on staggered grid, else = 0
    integer, dimension(:,:),pointer :: floating_mask => null()   !> = 1 where ice is present and floating, else = 0
    integer, dimension(:,:),pointer :: grounded_mask => null()   !> = 1 where ice is present and grounded, else = 0

    ! marine connection mask
    integer, dimension(:,:),pointer :: marine_connection_mask => null()   !> = 1 for cells with a marine connection to the ocean, else = 0

    ! masks for ice cap removal
    integer, dimension(:,:),pointer :: ice_sheet_mask => null()  !> = 1 for ice sheet cells, = 0 for ice cap cells
    integer, dimension(:,:),pointer :: ice_cap_mask => null()    !> = 1 for ice cap cells disconnected from the main ice sheet

    ! fields for forced retreat, as in ISMIP6 ocean-forced retreat for Greenland
    real(dp),dimension(:,:),pointer :: ice_fraction_retreat_mask => null() !> = 0.0 where ice is forced to retreat, = 1.0 where 
                                                                 !> ice is left intact; ice is thinned based on reference_thck for 
                                                                 !> values between 0 and 1
    real(dp),dimension(:,:),pointer :: reference_thck => null()  !> reference thickness giving upper limit for retreating ice

    ! fields for ice removal, e.g., removal of ice caps
    real(dp),dimension(:,:),pointer :: removal_thck => null()        !> thickness loss of ice removed in grid cell
    real(dp),dimension(:,:),pointer :: removal_rate => null()        !> rate of ice removal (m/yr ice)
    real(dp),dimension(:,:),pointer :: removal_rate_tavg => null()   !> rate of ice removal (m/yr ice, time average)

  end type glissade_geometry_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_scalars_type

    !> Holds various global and basin-scale diagnostic scalars that can be written to output

    ! global scalars
    real(dp) :: iarea                  ! total ice area (m^2)
    real(dp) :: iareag                 ! total grounded ice area (m^2)
    real(dp) :: iareaf                 ! total floating ice area (m^2)
    real(dp) :: ivol                   ! total ice volume (m^3)
    real(dp) :: ivol_above_flotation   ! total ice volume above flotation (m^3)
    real(dp) :: imass                  ! total ice mass (kg)
    real(dp) :: imass_above_flotation  ! total ice mass above flotation (kg)
    real(dp) :: icap_area              ! total ice cap area (m^2) (included within iarea)
    real(dp) :: icap_vol               ! total ice cap volume (m^3) (included within ivol)
    real(dp) :: total_bmlt_float          ! global sum of bmlt_float (kg s^-1)
    real(dp) :: total_bmlt_float_target   ! global sum of bmlt_float_target (kg s^-1)

    ! basin-scale scalars
    real(dp), dimension(:), pointer :: iarea_basin    ! total ice area per basin (m^2)
    real(dp), dimension(:), pointer :: iareag_basin   ! total grounded ice area per basin (m^2)
    real(dp), dimension(:), pointer :: iareaf_basin   ! total floating ice area per basin (m^2)
    real(dp), dimension(:), pointer :: ivol_basin     ! total ice volume per basin (m^3)
    real(dp), dimension(:), pointer :: ivol_above_flotation_basin   ! total ice volume above flotation per basin (m^3)
    real(dp), dimension(:), pointer :: imass_basin    ! total ice mass per basin (kg)
    real(dp), dimension(:), pointer :: imass_above_flotation_basin   ! total ice mass above flotation per basin (kg)
    real(dp), dimension(:), pointer :: icap_area_basin   ! total ice cap area per basin (m^2)
    real(dp), dimension(:), pointer :: icap_vol_basin    ! total ice cap volume per basin (m^3)
    real(dp), dimension(:), pointer :: bmlt_float_basin         ! total bmlt_float per basin (kg s^-1)
    real(dp), dimension(:), pointer :: bmlt_float_target_basin  ! total bmlt_float_target per basin (kg s^-1)

 end type glissade_scalars_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_mass_flux_type

    !> Holds various mass flux diagnostics

    ! mass fluxes at upper, lower and lateral boundaries
    ! Note: sfc_mbal_flux and basal_mbal_flux are equal to acab_applied and bmlt_applied multipled by a constant.
    ! Icebergs and isthmuses are included in the calving flux.
    ! The removal flux accounts for cells removed as ice caps without explicitly calving.

    real(dp),dimension(:,:),  pointer :: sfc_mbal_flux =>null()        !> surface mass balance (kg m^-2 s^-1), diagnosed from climate%acab_applied
    real(dp),dimension(:,:),  pointer :: sfc_mbal_flux_tavg =>null()   !> surface mass balance (kg m^-2 s^-1, time average)
    real(dp),dimension(:,:),  pointer :: basal_mbal_flux =>null()      !> basal mass balance (kg m^-2 s^-1), diagnosed from climate%bmlt_applied
    real(dp),dimension(:,:),  pointer :: basal_mbal_flux_tavg =>null() !> basal mass balance (kg m^-2 s^-1, time average)
    real(dp),dimension(:,:),  pointer :: calving_flux =>null()         !> calving flux (kg m^-2 s^-1), diagnosed from calving%calving_thck
    real(dp),dimension(:,:),  pointer :: calving_flux_tavg =>null()    !> calving flux (kg m^-2 s^-1, time average)
    real(dp),dimension(:,:),  pointer :: latmelt_flux =>null()         !> lateral melt flux (kg m^-2 s^-1), diagnosed from lateral_melt%melt_thck
    real(dp),dimension(:,:),  pointer :: latmelt_flux_tavg =>null()    !> lateral melt flux (kg m^-2 s^-1, time average)
    real(dp),dimension(:,:),  pointer :: removal_flux =>null()         !> removal flux (kg m^-2 s^-1), diagnosed from calving%%removal_thck
    real(dp),dimension(:,:),  pointer :: removal_flux_tavg =>null()    !> removal flux (kg m^-2 s^-1, time average)
    real(dp),dimension(:,:),  pointer :: gl_flux_east =>null()         !> mass flux eastward at grounding line, edge-based (kg m^-1 s^-1)
    real(dp),dimension(:,:),  pointer :: gl_flux_north =>null()        !> mass flux northward at grounding line, edge_based (kg m^-1 s^-1)
    real(dp),dimension(:,:),  pointer :: gl_flux =>null()              !> mass flux at grounding line, cell-based (kg m^-1 s^-1)
    real(dp),dimension(:,:),  pointer :: gl_flux_tavg =>null()         !> mass flux at grounding line, cell-based (kg m^-1 s^-1, time average)

    ! global scalars
    real(dp) :: total_smb_flux         ! total surface mass balance flux (kg/s)
    real(dp) :: total_smb_flux_tavg    ! total surface mass balance flux (kg/s), time average
    ! SMB budget diagnostics for coupled runs: total_smb_cpl_flux = total_smb_flux
    !  + total_smb_unapplied_glad_flux + total_smb_unapplied_flux
    real(dp) :: total_smb_cpl_flux = 0.d0          ! total SMB flux received from the coupler (kg/s)
    real(dp) :: total_smb_cpl_flux_tavg = 0.d0     ! total SMB flux received from the coupler (kg/s), time average
    real(dp) :: total_smb_unapplied_glad_flux = 0.d0 ! total SMB flux zeroed by GLAD (kg/s)
    real(dp) :: total_smb_unapplied_glad_flux_tavg = 0.d0 ! total SMB flux zeroed by GLAD (kg/s), time average
    real(dp) :: total_smb_unapplied_flux = 0.d0    ! total SMB flux passed to the dycore but not applied (kg/s)
    real(dp) :: total_smb_unapplied_flux_tavg = 0.d0 ! total SMB flux passed to the dycore but not applied (kg/s), time average
    real(dp) :: total_bmb_flux         ! total basal mass balance flux (kg/s)
    real(dp) :: total_bmb_flux_tavg    ! total basal mass balance flux (kg/s), time average
    real(dp) :: total_calving_flux     ! total calving mass flux (kg/s)
    real(dp) :: total_calving_flux_tavg! total calving mass flux (kg/s), time average
    real(dp) :: total_latmelt_flux     ! total latmelt mass flux (kg/s)
    real(dp) :: total_latmelt_flux_tavg! total latmelt mass flux (kg/s), time average
    real(dp) :: total_removal_flux     ! total removal mass flux (kg/s)
    real(dp) :: total_removal_flux_tavg! total removal mass flux (kg/s), time average
    real(dp) :: total_gl_flux          ! total grounding line mass flux (kg/s)
    real(dp) :: total_gl_flux_tavg     ! total grounding line mass flux (kg/s), time average

 end type glissade_mass_flux_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_geomderv_type

    !> Holds the horizontal and temporal derivatives of the thickness and
    !> upper surface elevation, as well as the thickness on the staggered grid.

    !*tb* Added a bunch of stuff here to clean up the higher order code that
    !I've been writing.  Might be worth it to add a mechanism to conditionally
    !allocate these depending on whether they are needed by the SIA core or by
    !the higher-order extensions

    !First derivatives on a staggered grid
    real(dp),dimension(:,:),pointer :: dthckdew => null() !> E-W derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdew => null() !> E-W derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdns => null() !> N-S derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdns => null() !> N-S derivative of upper surface elevation.

  end type glissade_geomderv_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_tensor_type
    real(dp), dimension(:,:,:), pointer :: scalar => null()
    real(dp), dimension(:,:,:), pointer :: xz => null()
    real(dp), dimension(:,:,:), pointer :: yz => null()
    real(dp), dimension(:,:,:), pointer :: xx => null()
    real(dp), dimension(:,:,:), pointer :: yy => null()
    real(dp), dimension(:,:,:), pointer :: xy => null()
  end type glissade_tensor_type
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_velocity_type

    !> Holds the velocity fields in 2D and 3D.
    real(dp),dimension(:,:,:),pointer :: uvel  => null()   !> 3D $x$-velocity.
    real(dp),dimension(:,:,:),pointer :: vvel  => null()   !> 3D $y$-velocity.
    real(dp),dimension(:,:,:),pointer :: velnorm => null() ! horizontal ice speed
    real(dp),dimension(:,:,:),pointer :: wvel  => null()   !> 3D $z$-velocity.

    ! Note: DIVA solves for uvel_2d and vvel_2d; these are typically (but not necessarily) the vertical average
    !       Velocities have units of m/s within the code.
    ! For I/O, velocities have units of m/yr, with a scale factor of 1/scyr
    real(dp),dimension(:,:)  ,pointer :: uvel_2d  => null()   !> 2D $x$-velocity (m/s); typically the vertical average
    real(dp),dimension(:,:)  ,pointer :: vvel_2d  => null()   !> 2D $y$-velocity (m/s); typically the vertical average
    real(dp),dimension(:,:)  ,pointer :: uvel_mean  => null() !> vertical mean $x$-velocity
    real(dp),dimension(:,:)  ,pointer :: vvel_mean  => null() !> vertical mean $y$-velocity 
    real(dp),dimension(:,:)  ,pointer :: velnorm_mean  => null() !> vertical mean ice speed

    real(dp),dimension(:,:)  ,pointer :: usfc_obs => null()     !> observed surface $x$-velocity
    real(dp),dimension(:,:)  ,pointer :: vsfc_obs => null()     !> observed surface $y$-velocity
    real(dp),dimension(:,:)  ,pointer :: velo_sfc_obs => null() !> observed surface speed = sqrt(usfc_obc^2 + vsfc_obs^2)
    real(dp),dimension(:,:)  ,pointer :: velo_sfc => null()     !> surface speed

    ! Note: uvel_extend and vvel_extend can be used for input and output of uvel, vvel on a staggered grid 
    !       that is the same size as the unstaggered grid. This is required for exact restart if velocities
    !       are nonzero along the north and east boundaries of the global domain.
    real(dp),dimension(:,:,:),pointer :: uvel_extend => null()  !> 3D $x$-velocity on extended staggered grid
    real(dp),dimension(:,:,:),pointer :: vvel_extend => null()  !> 3D $y$-velocity on extended staggered grid
    real(dp),dimension(:,:)  ,pointer :: uvel_2d_extend => null()  !> 2D $x$-velocity on extended staggered grid
    real(dp),dimension(:,:)  ,pointer :: vvel_2d_extend => null()  !> 2D $y$-velocity on extended staggered grid

    real(dp),dimension(:,:)  ,pointer :: btrc  => null()        !> basal traction (scalar field)
    real(dp),dimension(:,:)  ,pointer :: beta  => null()        !> basal shear coefficient on velo grid (Pa yr/m)
    real(dp),dimension(:,:)  ,pointer :: beta_internal => null()!> beta weighted by f_ground or otherwise adjusted (glissade only)
    real(dp),dimension(:,:)  ,pointer :: unstagbeta  => null()  !> basal shear coefficient on ice grid (Pa yr/m)

    !> mask that specifies where the velocity being read in should be held constant as a dirichlet condition
    integer, dimension(:,:), pointer  :: kinbcmask => null()    

    !> masks that specify where the outflow velocities on the global boundary should be set to zero
    integer, dimension(:,:), pointer  :: umask_no_penetration => null()
    integer, dimension(:,:), pointer  :: vmask_no_penetration => null()

    ! for viewing the spatial pattern of residuals
    real(dp),dimension(:,:,:),pointer :: resid_u => null()     ! u component of residual Ax - b where x is the velocity
    real(dp),dimension(:,:,:),pointer :: resid_v => null()     ! v component of residual Ax - b where x is the velocity

    ! for viewing the driving stress on the RHS
    real(dp),dimension(:,:,:),pointer :: rhs_u => null()       ! u component of b in Ax = b
    real(dp),dimension(:,:,:),pointer :: rhs_v => null()       ! v component of b in Ax = b

    ! strain rate tensor (HO only) 
    type(glissade_tensor_type) :: strain_rate                          ! strain rate tensor, diagnosed from stress tensor and efvs (s^-1)

    real(dp), dimension(:,:), pointer :: divu => null()        ! horizontal divergence rate, eps_xx + eps_yy (s^-1)
    real(dp), dimension(:,:), pointer :: shear => null()       ! invariant related to shear (s^-1): sqrt{[(eps_xx - eps_yy)/2]^2 + eps_xy^2}

    !> [velocity] keys

    integer :: stokes_approx = STOKES_APPROX_DIVA
    !> Flag that indicates which Stokes approximation to use with the glissade dycore.
    !> Not valid for other dycores 
    !> The default is diva; it was blatter_pattyn until the config reorganization.
    !> Note: There are two SIA options:
    !>       local_sia uses module glissade_velo_sia to compute local SIA velocities, similar to Glide
    !>       sia uses module glissade_velo_higher to compute SIA velocities via an iterative solve
    !> \begin{description}
    !> \item[local_sia] Shallow-ice approximation, Glide-type calculation; uses glissade_velo_sia
    !> \item[sia]  Shallow-ice approximation, vertical-shear stresses only; uses glissade_velo_higher
    !> \item[ssa]  Shallow-shelf approximation, horizontal-plane stresses only; uses glissade_velo_higher
    !> \item[blatter_pattyn]  Blatter-Pattyn approximation with both vertical-shear and horizontal-plane stresses; uses glissade_velo_higher
    !> \item[l1l2]  Vertically integrated 'L1L2' approximation with vertical-shear and horizontal-plane stresses; uses glissade_velo_higher
    !> \item[diva]  Depth-integrated viscosity approximation (DIVA) based on Goldberg (2011); uses glissade_velo_higher
    !> \item[hybrid]  Hybrid solver combining an SSA basal solve with a local vertical SIA solve
    !> \end{description}

    logical :: diva_slope_correction = .true.
    !> Flag that indicates whether the DIVA solver includes a correction for steep slopes

    integer :: efvs_method = EFVS_METHOD_NONLINEAR
    !> Flag that indicates how effective viscosity is computed
    !> \begin{description}
    !> \item[constant] constant value
    !> \item[flow_factor] multiple of flow factor
    !> \item[nonlinear] compute from effective strain rate
    !> \end{description}

    integer :: surface_gradient_method = SURFACE_GRADIENT_METHOD_CENTERED    
    !> Flag that indicates which gradient operator to use in the glissade dycore.
    !> Not valid for other dycores
    !> NOTE: Upstream may be better for ice evolution because it damps checkerboard noise.
    !> \begin{description}
    !> \item[centered] Centered gradient
    !> \item[upstream1] First-order accurate upstream gradient
    !> \item[upstream2] Second-order accurate upstream gradient
    !> \end{description}

    integer :: surface_gradient_margin_method = SURFACE_GRADIENT_MARGIN_METHOD_HYBRID
    !> Flag that indicates how to compute the gradient at the ice margin in the glissade dycore.
    !> Note: Gradients are always computed at edges with ice on both sides.
    !>       The methods differ in whether gradients are computed when an ice-covered cell
    !>        lies above an ice-free cell (land or ocean).
    !        Option 1 is appropriate for ice sheets with both land and marine boundaries,
    !         when a lateral spreading force is computed for marine ice cliffs.
    !         This lateral force ensures spreading where the gradient is zero.
    !> \begin{description}
    !> \item[land] Compute edge gradient when either cell is ice-covered
    !> \item[hybrid] Compute edge gradient for ice-covered cell above ice-free land (not ocean)
    !> \item[marine] Compute edge gradient only when both cells have ice
    !> \end{description}

    real(dp) :: efvs_const = 2336041.d0    ! value of efvs to use in constant efvs case, in units Pa yr
                                       ! = 0.5*A^(-1), where A = 2.140373 Pa^(-1) yr^(1) is the value used in ISMIP-HOM Test F
    real(dp) :: effstrain_min = 1.0d-8 ! min value of effective strain (yr^-1) for computation of effective viscosity
                                       ! 1.d-8 generally works well for SSA and DIVA, but a larger value (~1.e-6) may be needed for BP
    real(dp) :: max_slope = 1.0d0      ! maximum surface slope allowed in Glissade dycore (unitless)
                                       ! Note: It may be necessary to reduce max_slope to ~0.1 to prevent huge velocities
                                       !       in regions of rough coastal topography

  end type glissade_velocity_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_velocity_solver_type

    !> [velocity_solver] keys

    integer :: nonlinear_method = NONLINEAR_METHOD_PICARD
    !> Flag that indicates method for solving the nonlinear iteration when solving 
    !> the first-order momentum balance
    !> \begin{description}
    !> \item[picard] do standard Picard iteration
    !> \item[accelerated_picard] do Picard iteration with acceleration
    !> \end{description}

    integer :: nonlinear_residual = NONLINEAR_RESIDUAL_L2NORM
    !> Flag that indicates method for computing residual in PP dyn core: 
    !> \begin{description}
    !> \item[max] maxval 
    !> \item[max_no_ubas] maxval ignoring basal velocity 
    !> \item[mean] mean value
    !> \item[l2norm] L2 norm of system residual, Ax-b=resid
    !> \item[relative_l2norm] L2 norm of system residual relative to rhs, |Ax-b|/|b|
    !> \end{description}

    integer :: nonlinear_maxiters = 100    
    !> maximum number of nonlinear iterations to be used by the Glissade velocity solver

    integer :: linear_solver = LINEAR_SOLVER_PCG_CHRONGEAR
    !> Flag that indicates method for solving the sparse linear system
    !> that arises from the higher-order solver
    !> \begin{description}
    !> \item[slap_pcg] SLAP (serial): Preconditioned conjugate gradient, incomplete Cholesky preconditioner
    !> \item[slap_bicg]  SLAP (serial): Biconjugate gradient, incomplete LU preconditioner
    !> \item[slap_gmres]  SLAP (serial): GMRES, incomplete LU preconditioner
    !> \item[pcg] Native PCG, parallel-enabled, standard solver
    !> \item[pcg_chrongear] Native PCG, parallel-enabled, Chronopoulos-Gear solver (default)
    !> \end{description}

    integer :: preconditioner = PRECONDITIONER_DIAGONAL   ! default was PRECONDITIONER_SIA until the config reorganization
    !> Flag that indicates which Stokes preconditioner to use in the glissade dycore.
    !> Not valid for other dycores 
    !> \begin{description}
    !> \item[none] No preconditioner
    !> \item[diagonal] Diagonal preconditioner
    !> \item[sia] Physics-based shallow-ice preconditioner
    !> \item[tridiag_local] Local tridiagonal preconditioner (one solve per task)
    !> \item[tridiag_global] Global tridiagonal preconditioner (global solve)
    !> \end{description}

    integer :: linear_solve_ncheck = 5
    !> check the linear solver for convergence every linear_solve_ncheck iterations

    integer :: linear_maxiters = 200
    !> max number of linear iterations before quitting

    real(dp) :: linear_tolerance = 1.0d-08
    !> error tolerance for linear solver

    integer :: assemble_beta = ASSEMBLE_BETA_LOCAL
    !> Flag that describes how beta terms are assembled in the glissade finite-element calculation
    !> \begin{description}
    !> \item[standard] standard finite-element calculation (which effectively smooths beta at discontinuities)
    !> \item[local] apply local value of beta at each vertex
    !> \end{description}

    integer :: assemble_taud = ASSEMBLE_TAUD_LOCAL
    !> Flag that describes how driving-stress terms are assembled in the glissade finite-element calculation
    !> \begin{description}
    !> \item[standard] standard finite-element calculation (which effectively smooths the driving stress)
    !> \item[local] apply local value of driving stress at each vertex
    !> \end{description}

    integer :: assemble_bfric = ASSEMBLE_BFRIC_LOCAL
    !> Flag that describes how the basal friction heat flux is assembled in the glissade finite-element calculation
    !> \begin{description}
    !> \item[standard] standard finite-element calculation summing over quadrature points
    !> \item[local] apply local value of beta*(u^2 + v^2) at each vertex
    !> \end{description}

    integer :: assemble_lateral = ASSEMBLE_LATERAL_LOCAL   ! default was STANDARD until the config reorganization
    !> Flag that describes how lateral stress terms are assembled in the glissade finite-element calculation
    !> \begin{description}
    !> \item[standard] standard finite-element calculation of thck and usrf at quadrature points
    !> \item[local] apply local cell-center value of thck and usrf on each face
    !> \end{description}

  end type glissade_velocity_solver_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_advection_type

    !> [advection] keys

    integer :: advection_method = ADVECTION_METHOD_REMAP
    !> Thickness evolution method:
    !> \begin{description}
    !> \item[remap] Incremental remapping
    !> \item[upwind] 1st-order upwind scheme
    !> \item[fixed_thickness] Temperature evolves but thickness does not
    !> \end{description}

    integer :: vertical_remap_order = VERTICAL_REMAP_ORDER_FIRST
    !> Flag that indicates the order of accuracy for vertical remapping
    !> \begin{description}
    !> \item[first] first-order accurate in the vertical direction
    !> \item[second] second-order accurate in the vertical direction
    !> \end{description}

    logical :: compute_ice_age = .true.
    !> if true, then compute the ice age as a tracer

  end type glissade_advection_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_flow_law_type

    !> [flow_law] keys (n_glen is a module variable in glimmer_physcon)

    integer :: flow_law = FLOW_LAW_PATERSON_BUDD
    !> Method for calculating flow factor $A$:
    !> \begin{description} 
    !> \item[constant] Set equal to $1\times 10^{-16} \mathrm{Pa}^{-n} \mathrm{yr}^{-1}$
    !> \item[paterson_budd_minus5c] \emph{Paterson and Budd} relationship, 
    !> with temperature set to $-5^{\circ}\mathrm{C}$ 
    !> \item[paterson_budd] \emph{Paterson and Budd} relationship
    !> \end{description}

    logical :: damage_flwa_feedback = .false.
    !> if true, then flwa increases as damage increases
    !> This can be a strong positive feedback on damage

    real(dp) :: flow_factor_ground = 1.0d0 ! flow enhancement parameter for the Arrhenius relationship;
                                           ! grounded ice only; typically > 1 for SIA models to speed up the ice
    real(dp) :: flow_factor_float = 1.0d0  ! flow enhancement parameter for floating ice
                                           ! Default is 1.0, but for marine simulations a smaller value
                                           !  is often needed to match observed shelf speeds
    real(dp) :: flwa_const = 1.0d-16       ! Glen's A to use in isothermal case, in units Pa^{-n} yr^{-1} 
                                           ! (would change to e.g. 4.6e-18 in EISMINT-ROSS case)

  end type glissade_flow_law_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_grounding_line_type

    !> [grounding_line] keys

    integer :: grounding_line_param = GROUNDING_LINE_PARAM_NONE
    !> Flag that indicates how to compute the grounded fraction of each gridcell in the glissade dycore.
    !> Not valid for other dycores
    !> \begin{description}
    !> \item[none] f_ground = 0 in floating cells (based on flotation condition), else f_ground = 1
    !> \item[vertices] basal friction GLP with 0 <= f_ground <= 1 at vertices
    !> \item[deluxe] "deluxe" GLP with f_ground at both cells and vertices; includes basal fraction and basal melting
    !> \end{description}

    integer :: fground_no_glp = FGROUND_NO_GLP_GROUNDED_NEIGHBOR
    !> Flag that indicates how to identify grounded and floating vertices when running without a GLP
    !> \begin{description}
    !> \item[grounded_neighbor] f_ground = 1 if any neighbor cell is grounded or land, else f_ground = 0
    !> \item[flotation_function] f_ground = 0 or 1 based on staggered flotation function
    !> \end{description}

    integer :: grounding_line_bmlt = GROUNDING_LINE_BMLT_NONE
    !> Flag that indicates how to compute bmlt_float in partly grounded cells
    !> \begin{description}
    !> \item[none] Apply bmlt_float in all floating cells, based on floating_mask
    !> \item[floating_fraction] Weigh bmlt by floating fraction of cell if 0 < f_ground < 1
    !> \item[zero_grounded] Set bmlt_float = 0 in partly grounded cells (f_ground > 0)
    !> \end{description}

    integer :: flotation_function = FLOTATION_FUNCTION_LINEAR
    !WHL TODO - Remove option 3 during the config reorg
    !> Flag that indicates how to compute the flotation function at and near vertices in the glissade dycore
    !> Not valid for other dycores
    !> \begin{description}
    !> \item[pattyn] f_flotation = (-rhow*b/rhoi*H) = f_pattyn; <=1 for grounded, > 1 for floating
    !> \item[inverse_pattyn] f_flotation = (rhoi*H)/(-rhow*b) = 1/f_pattyn; >=1 for grounded, < 1 for floating
    !> \item[linear] f_flotation = -rhow*b - rhoi*H = ocean cavity thickness; <=0 for grounded, > 0 for floating 
    !>          (without extrapolation to ice-free ocean, as in an older version)
    !> \item[3] Like (2), but with a correction for topg
    !> \end{description}

  end type glissade_grounding_line_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_removal_type

    !> [removal] keys

    logical :: remove_ice_caps = .false.
    !> Flag that indicates whether ice caps are removed and added to the calving flux

    integer :: force_retreat = FORCE_RETREAT_NONE
    !> Flag that indicates whether retreat is forced using ice_fraction_retreat_mask
    !> item[0] do not force retreat
    !> item[1] force retreat of all ice identified by a retreat mask
    !> item[2] force retreat of floating or weakly grounded ice identified by a retreat mask

  end type glissade_removal_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_geothermal_type

    !> [geothermal] keys

    integer :: geothermal_method = GEOTHERMAL_METHOD_UNIFORM
    !> geothermal heat flux:
    !> \begin{description}
    !> \item[uniform] prescribed uniform geothermal heat flux
    !> \item[file] read 2D geothermal flux field from input file (if present)
    !> \end{description}

    real(dp) :: geothermal_flux = -5.0d-2  ! uniform geothermal flux (W m^{-2}), positive down  !TODO - Change default to 0.0?

  end type glissade_geothermal_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_stress_type      

    type(glissade_tensor_type) :: tau         ! HO only
    real(dp),dimension(:,:,:),pointer :: efvs => null()    !> effective viscosity
    real(dp),dimension(:,:),  pointer :: efvs_vertavg => null()  !> vertical mean effective viscosity
    real(dp),dimension(:,:),  pointer :: btractx => null() !> basal traction (Pa), x comp
    real(dp),dimension(:,:),  pointer :: btracty => null() !> basal traction (Pa), y comp
    real(dp),dimension(:,:),  pointer :: btract => null()  !> basal traction (Pa), magnitude = sqrt(btractx^2 + btracty^2)
    !WHL - The extended versions are needed for exact restart if using DIVA solver for a problem with nonzero traction at global boundaries
    real(dp),dimension(:,:),  pointer :: btractx_extend => null() !> basal traction (Pa), x comp, on extended staggered grid
    real(dp),dimension(:,:),  pointer :: btracty_extend => null() !> basal traction (Pa), y comp, on extended staggered grid
    real(dp),dimension(:,:),  pointer :: taudx => null()   !> driving stress (Pa), x comp
    real(dp),dimension(:,:),  pointer :: taudy => null()   !> driving stress (Pa), y comp

  end type glissade_stress_type      

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Make eus a config file parameter to support runs with time-varying sea level forcing

  ! Note: smb and acab are both SMB variables, but they have different units.
  !       smb has units of kg/m^2/yr w.e. = mm/yr w.e.
  !       acab has units of m/yr ice for I/O, and SI units in most of the code.
  !       Terms with a prefix 'smb_' or 'acab_' follow the same units convention.

  ! Note on acab_applied_tavg: This is the average value of acab_applied over an output interval.
  !                    If 'average = 1' in the acab_applied entry of glissade_vars.def, then acab_applied_tavg is automatically
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

  type glissade_climate_type
     !> Holds fields used to drive the model
     real(dp),dimension(:,:),pointer :: acab            => null() !> Surface mass balance (m/yr ice for I/O, m/s ice in code)
     real(dp),dimension(:,:),pointer :: acab_applied    => null() !> Surface mass balance applied to ice (m/yr ice)
                                                                  !>    = 0 for ice-free cells with acab < 0
     real(dp),dimension(:,:),pointer :: acab_applied_tavg => null() !> Surface mass balance applied to ice (m/yr ice, time average)
     ! SMB budget diagnostics for coupled runs (GLAD); m/yr ice for I/O, m/s ice in code
     ! In each cell, acab_cpl = acab_applied + acab_unapplied_glad + acab_unapplied
     real(dp),dimension(:,:),pointer :: acab_cpl        => null() !> SMB received from the coupler, before any masking in GLAD
     real(dp),dimension(:,:),pointer :: acab_cpl_tavg   => null() !> SMB received from the coupler (time average)
     real(dp),dimension(:,:),pointer :: acab_unapplied_glad      => null() !> SMB zeroed by GLAD before the dycore sees it
                                                                           !>    (ice-free cells with acab < 0; ice-free ocean)
     real(dp),dimension(:,:),pointer :: acab_unapplied_glad_tavg => null() !> SMB zeroed by GLAD (time average)
     real(dp),dimension(:,:),pointer :: acab_unapplied      => null() !> SMB passed to the dycore but not applied
                                                                      !>    (melt exceeding the ice present; ocean cells; partial cells;
                                                                      !>     positive SMB removed by block_inception)
     real(dp),dimension(:,:),pointer :: acab_unapplied_tavg => null() !> SMB passed to the dycore but not applied (time average)
     real(dp),dimension(:,:),pointer :: smb             => null() !> Surface mass balance (mm/yr water equivalent)
     real(dp),dimension(:,:),pointer :: smb_anomaly     => null() !> Surface mass balance anomaly (mm/yr w.e.)
     real(dp),dimension(:,:),pointer :: smb_corrected  => null()  !> Surface mass balance with flux or anomaly corrections (mm/yr w.e.)
     real(dp),dimension(:,:),pointer :: artm            => null() !> Annual mean air temperature (degC)
     real(dp),dimension(:,:),pointer :: artm_anomaly    => null() !> Annual mean air temperature anomaly (degC)
     real(dp),dimension(:,:),pointer :: artm_corrected  => null() !> Annual mean air temperature with anomaly corrections (degC)
     real(dp),dimension(:,:),pointer :: snow            => null() !> snowfall rate (mm/yr w.e.)
     real(dp),dimension(:,:),pointer :: snow_anomaly    => null() !> snowfall anomaly (mm/yr w.e.)
     real(dp),dimension(:,:),pointer :: snow_corrected  => null() !> snowfall with anomaly corrections (mm/yr w.e.)
     real(dp),dimension(:,:),pointer :: precip          => null() !> precipitation rate (mm/yr w.e.)
                                                                  !> for glaciers, snow can be derived from precip + downscaled artm
     real(dp),dimension(:,:),pointer :: precip_anomaly  => null() !> precip anomaly (mm/yr w.e.)
     real(dp),dimension(:,:),pointer :: precip_corrected=> null() !> precip with anomaly corrections (mm/yr w.e.)
     integer, dimension(:,:),pointer :: overwrite_acab_mask => null() !> mask for cells where acab is overwritten
     real(dp),dimension(:,:),pointer :: smb_obs         => null() !> Observed surface mass balance (mm/yr water equivalent)
                                                                  !> 'smb' could have any source (models, obs, etc.), but smb_obs
                                                                  !>  is always from observations and may be an inversion target

     ! Next several fields used for SMB_INPUT_FUNCTION_GRADZ, ARTM_INPUT_FUNCTION_GRADZ, ARTM_INPUT_FUNCTION_LAPSE
     ! Note: If both smb and artm are input in this format, they share the array usrf_ref.
     !       Sign convention for gradz is positive up, so artm_gradz is usually negative.
     real(dp),dimension(:,:),pointer :: smb_ref  => null()        !> SMB at reference elevation (mm/yr w.e.)
     real(dp),dimension(:,:),pointer :: smb_gradz  => null()      !> vertical gradient of SMB (mm/yr w.e. per m), positive up
     real(dp),dimension(:,:),pointer :: artm_ref => null()        !> artm at reference elevation (deg C)
     real(dp),dimension(:,:),pointer :: artm_gradz => null()      !> vertical gradient of artm (deg C per m), positive up
     real(dp),dimension(:,:),pointer :: usrf_ref => null()        !> reference upper surface elevation before lapse rate correction (m)

     ! Next several fields used for SMB_INPUT_FUNCTION_XYZ, ARTM_INPUT_FUNCTION_XYZ
     ! Note: If both smb and artm are input in this format, they share the array zatm(nzatm).
     integer  :: nzatm = 1                          !> number of atmosphere levels at which smb and other fields are provided
!     real(dp) :: dzatm = 0.d0                       !> could be used for evenly spaced layers; not currently used
     real(dp), dimension(:), pointer :: &
          zatm => null()                            !> atmosphere levels (m) where forcing is provided; positive up

     !Note: The first index is the vertical dimension; the second and third are the horizontal dimensions
     real(dp),dimension(:,:,:),pointer :: smb_3d        => null() !> SMB at multiple vertical levels (mm/yr w.e.)
     real(dp),dimension(:,:,:),pointer :: artm_3d       => null() !> artm at multiple vertical levels (m/yr ice)

     ! Ablation is computed independent of accumulation for SMB_INPUT_FUNCTION_PDD
     real(dp),dimension(:,:),pointer :: ablation => null()        !> surface ablation (mm/yr w.e.)

     ! The next several fields are used for the 'read_once' forcing option.
     ! E.g., if we want to read in all time slices of precip at once, we would set 'read_once' = .true. in the config file.
     ! All time slices are then stored in the precip_read_once array, where the third dimension is the number of time slices.
     ! Data are copied from precip_read_once to the regular 2D precip array as the model runs forward in time.
     real(dp), dimension(:,:,:), pointer :: &
          artm_ref_read_once => null(),           & !> artm_ref field, read_once version
          snow_read_once => null(),               & !> snow field, read_once version
          precip_read_once => null()                !> precip field, read_once version

     real(dp), dimension(:,:,:), pointer :: &
          artm_anomaly_read_once => null(),       & !> anomaly artm_ref field, read_once version
          snow_anomaly_read_once => null(),       & !> anomaly snow field, read_once version
          precip_anomaly_read_once => null()        !> anomaly precip field, read_once version

     real(dp) :: eus = 0.d0                         !> eustatic sea level

     real(dp) :: smb_anomaly_tstart = 0.0d0         !> time to start applying the anomaly (yr)
     real(dp) :: smb_anomaly_timescale = 0.0d0      !> number of years over which the smb anomaly is phased in linearly
                                                    !> If set to zero, then the anomaly is applied immediately.
                                                    !> The initMIP value is 40 yr.
     real(dp) :: overwrite_acab_value = 0.0d0       !> acab value to apply in grid cells where overwrite_acab_mask = 1
                                                    !> input as m/yr ice, converted to model units
     real(dp) :: overwrite_acab_thck_min = 0.0d0     !> overwrite acab where thck <= overwrite_acab_minthck
     real(dp) :: artm_anomaly_const = 0.0d0         !> spatially uniform value of artm_anomaly (degC)
     real(dp) :: artm_anomaly_tstart = 0.0d0        !> time to start applying the anomaly (yr)
     real(dp) :: artm_anomaly_timescale = 0.0d0     !> number of years over which the artm anomaly is phased in linearly
                                                    !> If set to zero, then the anomaly is applied immediately.
                                                    !> Snow and precip anomalies are assumed to have the same timescale
     real(dp) :: t_lapse = 0.0d0                    !> air temp lapse rate (deg/m); positive for T decreasing with height
     real(dp) :: degree_factor = 1000.d0            !> ablation rate as a function of temperature above tmlt (mm/yr w.e./deg C)
                                                    !> Divide by 365 to convert to mm/day
     real(dp) :: tmlt = -1.d0                       !> spatially uniform temperature threshold for melting (deg C)
     real(dp) :: snow_threshold_min = 0.0d0         !> air temperature (deg C) below which all precip falls as snow
     real(dp) :: snow_threshold_max = 2.0d0         !> air temperature (deg C) above which all precip falls as rain

     !TODO - Remove the equivalent glacier variables


    !> [climate] options (parameters are above)

    integer :: smb_input_units = SMB_INPUT_UNITS_M_ICE_PER_YR
    !> units for SMB input:
    !> \begin{description}
    !> \item[m_ice_per_yr] SMB input in units of m/yr ice (same as acab)
    !> \item[mm_we_per_yr] SMB input in units of mm/yr water equivalent
    !> \end{description}

    integer :: smb_input_function = SMB_INPUT_FUNCTION_XY
    !> functional form of SMB input:
    !> \begin{description}
    !> \item[xy] SMB(x,y); input as a function of horizontal location only
    !> \item[xy_gradz] SMB(x,y) + dSMB/dz(x,y) * dz; input SMB and its vertical gradient
    !> \item[xyz] SMB(x,y,z); input SMB at multiple elevations
    !> \item[pdd] SMB computed with a positive-degree scheme
    !> \end{description}

    integer :: artm_input_function = ARTM_INPUT_FUNCTION_XY
    !> functional form of surface temperature (artm) input:
    !> \begin{description}
    !> \item[xy] artm(x,y); input as a function of horizontal location only
    !> \item[xy_gradz] artm(x,y) + dartm/dz(x,y) * dz; input artm and its vertical gradient
    !> \item[xyz] artm(x,y,z); input artm at multiple elevations
    !> \item[xy_lapse] artm(x,y) + tlapse * dz; input artm and uniform lapse rate
    !> \end{description}

    logical :: enable_smb_anomaly = .false.
    !> if true, then apply a prescribed anomaly to smb

    logical :: enable_artm_anomaly = .false.
    !> if true, then apply a prescribed anomaly to artm

    logical :: enable_snow_anomaly = .false.
    !> if true, then apply a prescribed anomaly to snow

    logical :: enable_precip_anomaly = .false.
    !> if true, then apply a prescribed anomaly to precip

    integer :: overwrite_acab = OVERWRITE_ACAB_NONE
    !> overwrite acab (m/yr ice) in selected regions:
    !> \begin{description}
    !> \item[none] Do not overwrite acab anywhere
    !> \item[acab_zero] Overwrite acab where input acab = 0
    !> \item[thck_min] Overwrite acab where input thickness <= threshold value
    !> \item[mask] Overwrite acab where input mask = 1
    !> \end{description}

    logical :: enable_acab_dthck_dt_correction = .false.
    !> if true, then add (-dthck_dt_obs) to acab for floating ice

    logical :: block_inception = .false.
    !> Flag that indicates whether ice inception away from the main ice sheet is blocked

  end type glissade_climate_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_calving_type
     !> holds fields and parameters related to calving
     real(dp),dimension(:,:),  pointer :: calving_thck => null()   !> thickness loss in grid cell due to calving during one time step (m)
     real(dp),dimension(:,:),  pointer :: calving_rate => null()   !> rate of ice loss due to calving (m/yr ice)
     real(dp),dimension(:,:),  pointer :: calving_rate_tavg => null()  !> rate of ice loss due to calving (m/yr ice, time average)
     integer, dimension(:,:),  pointer :: calving_mask => null()         !> calve floating ice where the mask = 1 (calving_law = CALVING_LAW_GRID_MASK)
     real(dp),dimension(:,:),  pointer :: subgrid_calving_mask => null() !> calve floating ice where the mask < 1.0 (calving_law = CALVING_LAW_GRID_MASK);
                                                                         !> real instead of integer for use with subgrid calving parameterization
     integer, dimension(:,:),  pointer :: calving_front_mask => null()   !> = 1 for cells on the calving front, else = 0
     integer, dimension(:,:),  pointer :: beyond_cf_mask => null() !> = 1 for cells beyond the CF when using the subgrid CF scheme;
                                                                   !> these cells not allowed to fill until upstream neighbors are full
     real(dp),dimension(:,:),  pointer :: thck_effective => null() !> effective thickness for calving (m)
     real(dp),dimension(:,:),  pointer :: effective_areafrac => null() !> effective fractional area, < 1 for partial CF cells (m)
     real(dp),dimension(:,:),  pointer :: tau_eigen1 => null()     !> first eigenvalue of 2D horizontal stress tensor (Pa)
     real(dp),dimension(:,:),  pointer :: tau_eigen2 => null()     !> second eigenvalue of 2D horizontal stress tensor (Pa)
     real(dp),dimension(:,:),  pointer :: eps_eigen1 => null()     !> first eigenvalue of 2D horizontal strain rate tensor (s^-1)
     real(dp),dimension(:,:),  pointer :: eps_eigen2 => null()     !> second eigenvalue of 2D horizontal strain rate tensor (s^-1)
     real(dp),dimension(:,:,:),pointer :: damage => null()         !> 3D damage tracer, 0 > damage < 1 (calving_law = CALVING_LAW_DAMAGE)
  
     real(dp) :: marine_limit =  -200.d0         !> value of topg/relx at which floating ice calves (m)
                                                 !> (calving_law = CALVING_LAW_RELX_THRESHOLD, CALVING_LAW_TOPG_THRESHOLD)
     real(dp) :: calving_timescale = 0.0d0               !> timescale (yr) for calving (Glissade only); calving_thck = thck*min(dt/calving_timescale,1)
                                                 !> if calving_timescale = 0, then the full column calves at once
     real(dp) :: calving_thck_threshold = 0.d0                  !> minimum thickness (m) of floating ice at marine edge before it calves;
                                                 !> if used, must be set to a nonzero value in the config file
                                                 !> (calving_law = CALVING_LAW_THCK_THRESHOLD, EIGENCALVING, CALVING_LAW_DAMAGE)
     real(dp) :: dthck_dx_cf = 0.0d0             !> assumed max value of |dH/dx| at the calving front for full (not partial) cells (m/m)
     real(dp) :: thck_effective_min = 50.0d0     !> minimum value of thck_effective (m) for calving-front cells
     real(dp) :: eigencalving_const = 0.0d0           !> constant that multiplies the eigen-based strain rate to determine the calving rate (m)
     real(dp) :: tau_eigenconstant1 = 0.25d0     !> constant that multiplies the tau_eigen1 term in stress-based calving (unitless)
     real(dp) :: tau_eigenconstant2 = 1.0d0      !> constant that multiplies the tau_eigen2 term in stress-based calving (unitless)
     real(dp) :: stress_threshold = 6.0e4        !> stress threshold for CF retreat with stress-based calving (Pa)
     real(dp) :: lateral_rate_min = 0.0d0        !> minimum lateral retreat rate at CF (m/yr)
     real(dp) :: effec_stress_min = 0.0d0        !> minimum effective stress for calving at CF (Pa)
     real(dp) :: calving_length_scale = 4000.d0          !> length scale (m) for calving events (CALVING_LAW_STRESS_STOCHASTIC)
     real(dp) :: damage_const = 1.e5          !> damage scaling term with units of Pa * yr
!     real(dp) :: damage_constant1 = 0.0d0        !> damage constant that multiplies tau_eigen1 (yr^-1)
!     real(dp) :: damage_constant2 = 0.0d0        !> damage constant that multiplies tau_eigen2 (yr^-1)
     real(dp) :: taumax_cliff = 1.0d6            !> yield stress (Pa) for marine-based ice cliffs
     real(dp) :: cliff_timescale = 10.0d0        !> time scale (yr) for limiting marine cliffs (yr)
     real(dp) :: calving_front_x = 0.0d0         !> for options with a calving mask, calve ice wherever abs(x) > calving_front_x (m)
     real(dp) :: calving_front_y = 0.0d0         !> for options with a calving mask, calve ice wherever abs(y) > calving_front_y (m)
     real(dp) :: calving_front_radius = 0.0d0    !> for options with a calving mask, calve ice where the distance from the origin > radius
                                                 !> NOTE: Applied only if calving_front_x, calving_front_y, or calving_front_radius > 0
     real(dp) :: f_ground_threshold = 0.10d0     !> Threshold fraction for grounded cells in iceberg removal algorithm
                                                 !> Also used for isthmus removal

     ! calvingMIP parameters and diagnostics
     ! Note: For the circular domain, axis 1 is the y-axis and axis 2 is the line y = x in the NE quadrant
     !       For the Thule domain, axis 1 is the Caprona A axis, and axis 2 is the Halbrane A axis, both in the NW quadrant
     real(dp) :: &
          cf_advance_retreat_amplitude = 0.0d0,& !> prescribed amplitude (m/yr) for calving front advance or retreat
                                                 !> positive for sin(2*pi*t/period), negative for -sin(2*pi*t/period)
                                                 !> should be negative for CalvingMIP Experiments 2 and 4
          cf_advance_retreat_period = 0.0d0      !> period (yr) for an advance/retreat cycle
                                                 !> period = 0 => constant amplitude

     ! The following are for calvingMIP diagnostics along 8 axes
     ! Could be generalized for other problems with idealized geometry

     integer :: naxis = 1                                    !> number of axes for calvingMIP diagnostics
                                                             !> set to 8 below if running on a calvingMIP domain
     integer, dimension(:), pointer :: axis => null()        !> axis dimension variable, used for I/O

     real(dp), dimension(:), pointer :: cf_locx => null()    !> CF location, x coordinate (m) along each axis
     real(dp), dimension(:), pointer :: cf_locy => null()    !> CF location, y coordinate (m) along each axis
     real(dp), dimension(:), pointer :: cf_radius => null()  !> distance of calving front from origin (m) along each axis
     real(dp), dimension(:), pointer :: cf_thck => null()    !> ice thickness at CF (m) along each axis
     real(dp), dimension(:), pointer :: cf_uvel => null()    !> ice speed at CF (m/s), u component along each axis
     real(dp), dimension(:), pointer :: cf_vvel => null()    !> ice speed at CF (m/s), v component along each axis


    !> [calving] options (parameters are above)

    integer :: calving_law = CALVING_LAW_FLOAT_ZERO
    !> Calving: 
    !> \begin{description} 
    !> \item[none] No calving 
    !> \item[float_zero] Set thickness to zero if floating
    !> \item[relx_threshold] Set thickness to zero if relaxed bedrock is more than a
    !>          certain water depth (variable "marine_limit" in glissade_types)  
    !> \item[topg_threshold] Set thickness to zero if present bedrock topography lies below
    !>          a certain water depth (variable "marine_limit" in glissade_types)  
    !> \item[grid_mask] Calve based on grid location using a prescribed mask
    !>          (field 'calving_mask' or 'subgrid_calving_mask')
    !> \item[advance_retreat_rate] Prescribe the rate of calving front advance or retreat
    !> \item[thck_threshold] Calve ice whose thickness is below a given threshold
    !> \item[stress] Deterministic calving based on eigenvalues of the horizontal stress tensor
    !> \item[stress_stochastic] Stochastic calving based on eigenvalues of the horizontal stress tensor
    !> \item[damage] Calving based on stochastic internal damage
    !> \item[eigencalving] Eigencalving following Levermann et al. (2012)
    !> \end{description}

    logical :: calve_at_startup = .false.
    !> if true, then calve at initialization

    integer :: calving_domain = CALVING_DOMAIN_EVERYWHERE
    !> \begin{description}
    !> \item[ocean_edge] Calve only at ocean edge
    !> \item[everywhere] Calve wherever the calving criterion is met
    !> \end{description}

    logical :: apply_calving_mask = .false.
    !> if true, then apply a calving mask to prevent calving-front advance

    logical :: remove_icebergs = .true.
    !> if true, then identify and remove icebergs after calving
    !> These are connected regions with zero basal traction and no connection to grounded ice.
    !>       Safer to make it true, but not necessary for all applications

    logical :: remove_isthmuses = .false.
    !> if true, then identify and remove ice isthmuses after calving
    !> These are narrow bridges connecting two regions of floating ice.
    !> False by default, but may need to be true for the FORCE_RETREAT_FLOATING_ICE option.
    !> Requires remove_icebergs = .true. (for stability).

    logical :: limit_marine_cliffs = .false.
    !> if true, then thin marine-based cliffs based on a thickness threshold

    logical :: subgrid_calving_front = .false.
    !> if true, then use a subgrid calving front parameterization with partial cells (floating only) at the calving front

    integer :: calvingmip_domain = CALVINGMIP_DOMAIN_NONE
    !> Flag that indicates the desired domain for CalvingMIP experiments
    !> \begin{description}
    !> \item[none] none
    !> \item[circular] circular (radially symmetric)
    !> \item[circular] Thule (complex topography)
    !> \end{description}

  end type glissade_calving_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_lateral_melt_type

     ! holds parameters for lateral melting at floating or marine-grounded ice fronts

     real(dp),dimension(:,:),  pointer :: latmelt_thck => null()         !> thickness loss in grid cell due to lateral melt
     real(dp),dimension(:,:),  pointer :: latmelt_rate => null()         !> rate of ice loss due to lateral melt (m/yr ice)
     real(dp),dimension(:,:),  pointer :: latmelt_rate_tavg => null()    !> rate of ice loss due to lateral melt (m/yr ice, time average)

     integer ,dimension(:,:),  pointer :: melt_front_mask => null()   !> = 1 for cells on the lateral melt front, else = 0

     real(dp) :: melt_rate_const = 0.0d0   !> constant lateral retreat rate at melt front (m/yr)
     real(dp) :: melt_factor = 1.0d0       !> multiplier for Rignot frontal melt. A value of 1.6 was proposed for ISMIP7

     ! Can either read runoff_applied from an input file or compute it from surface ablation
     real(dp),dimension(:,:), pointer :: runoff_applied => null() !> runoff applied to subglacial discharge(kg/m^2/s)
     logical :: runoff_from_ablation = .false.  !> if true, then compute runoff_applied from surface ablation;
                                                !> if false, then read runoff_appied from an input file

     ! Can either read thermal_forcing_applied from an input file or compute it from 3d thermal forcing
     real(dp), dimension(:,:), pointer :: &
          thermal_forcing_applied => null()             !> thermal forcing applied to lateral melt, averaged over a depth range (deg K)
     logical :: thermal_forcing_avg_3d_to_2d = .false.  !> if true, then read 3d thermal forcing from file and average to 2d;
                                                        !> if false, then read thermal_forcing_applied from an input file
     real(dp) :: ztop_tfavg = -200.d0                  !> top end of depth range (m) for thermal_forcing_applied
     real(dp) :: zbot_tfavg = -500.d0                  !> bottom end of depth range (m) for thermal_forcing_applied


    !> [lateral_melt] switch

    integer :: lateral_melt_method = LATERAL_MELT_METHOD_NONE
    !> Lateral melt:
    !> \begin{description}
    !> \item[none] No lateral melt
    !> \item[constant] Constant lateral melt rate
    !> \item[ismip] ISMIP lateral melt rate based on thermal forcing and subglacial discharge

  end type glissade_lateral_melt_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_temperature_type

    !> Holds fields relating to temperature.

    !Note: With remapping advection of temperature, temp, flwa and dissip live on the
    !       staggered vertical grid at layer midpoints.
    !       The vertical dimensions are (0:upn) for temp and (1:upn-1) for flwa and dissip.
    !
    !      bheatflx and ucondflx are defined as positive down,
    !       so they will often be < 0.  
    !      However, bfricflx and bhydroflx are defined to be >= 0.
    !
    !      If bheatflx is read from a data file, be careful about the sign!
    !      In input data, the geothermal heat flux may be defined as positive upward,
    !       whereas bheatflx is defined as positive downward.

    real(dp),dimension(:,:,:),pointer :: temp => null()      !> 3D temperature field.
    real(dp),dimension(:,:,:),pointer :: tempunstag => null()!> 3D temperature field unstaggered in the vertical.
    real(dp),dimension(:,:),  pointer :: bheatflx => null()  !> basal heat flux (W/m^2) (geothermal, positive down)
    real(dp),dimension(:,:,:),pointer :: flwa => null()      !> Glen's flow factor $A$.
    real(dp),dimension(:,:,:),pointer :: dissip => null()    !> interior heat dissipation rate, divided by rhoi*Ci (deg/s)
    real(dp),dimension(:,:),  pointer :: btemp => null()     !> Basal temperature on ice grid; diagnosed from temp(upn)
    ! Note: With a basal temperature GLP, a grid cell can have different basal temps for grounded and floating ice
    real(dp),dimension(:,:),  pointer :: btemp_ground => null() !> Basal temperature for grounded ice
    real(dp),dimension(:,:),  pointer :: btemp_float => null()  !> Basal temperature for floating ice
    real(dp),dimension(:,:),  pointer :: bpmp => null()      !> Basal pressure melting point temperature
    real(dp),dimension(:,:),  pointer :: bfricflx => null()  !> basal heat flux (W/m^2) from friction (>= 0)
    real(dp),dimension(:,:),  pointer :: bhydroflx => null() !> basal heat flux (W/m^2) from refreezing meltwater (>= 0)
    real(dp),dimension(:,:,:),pointer :: waterfrac => null() !> fractional water content in layer (0 <= waterfrac <= 1)
    real(dp),dimension(:,:,:),pointer :: enthalpy => null()  !> specific enthalpy in layer (J m-3)
                                                             !> = rhoi * Ci * T for cold ice
    real(dp),dimension(:,:),  pointer :: ucondflx => null()  !> conductive heat flux (W/m^2) at upper sfc (positive down)

    real(dp),dimension(:,:),  pointer :: flow_enhancement_factor => null()  !> flow enhancement factor E (unitless)

    real(dp) :: pmp_offset = 5.0d0        ! offset of initial Tbed from pressure melting point temperature (deg C)
    real(dp) :: pmp_threshold = 1.0d-3    ! bed is assumed thawed where Tbed >= pmptemp - pmp_threshold (deg C)


    !> [temperature] keys

    integer :: temperature_method = TEMPERATURE_METHOD_PROGNOSTIC
    !> Method of ice temperature calculation:
    !> \begin{description} 
    !> \item[isothermal_artm] Set column to surface air temperature
    !> \item[prognostic] Prognostic temperature solution 
    !> \item[fixed] Do NOTHING - hold temperatures steady at initial value  
    !> \item[enthalpy] Prognostic enthalpy solution
    !> \end{description}

    integer :: temperature_init = TEMPERATURE_INIT_ARTM
    !> Temperature initialization:
    !> \begin{description} 
    !> \item[zero_celsius] Initialize temperature to 0 C
    !> \item[artm] Initialize temperature to surface air temperature
    !> \item[linear] Initialize temperature with a linear profile in each column
    !> \item[advective_diffusive] Initialize temperature with an advective-diffusive balance in each column
    !> \item[file] Initialize temperature from external file
    !> \end{description}

    integer :: dissipation_method = DISSIPATION_METHOD_FIRST_ORDER
    !> Flag that indicates method for computing the dissipation during the temperature calc.
    !> \begin{description}
    !> \item[none] for no dissipation
    !> \item[sia] for 0-order SIA approx
    !> \item[first_order] for first-order dissipation (Blatter-Pattyn)
    !> \end{description}

  end type glissade_temperature_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_inversion_type

     !TODO - Break into different derived types for each kind of inversion?

     ! parameters for initializing inversion fields
     real(dp) :: &
          inversion_thck_threshold = 0.0d0,          & !> ice thinner than this threshold (m), or than thcklim if larger, is removed at initialization
          inversion_thck_flotation_buffer = 1.0d0      !> if usrf_obs implies thck near the flotation thickness,
                                             !> set to thck_flotation +/- thck_flotation_buffer (m)
     ! general inversion parameters
     real(dp) :: &
          inversion_damping_factor = 2.0d0     !> factor that multiplies the dH/dt term in the inversion

     ! fields and parameters for powerlaw_c and coulomb_c inversion
     ! Note: powerlaw_c and coulomb_c are in the basal_physics type

     ! parameters for adjusting powerlaw_c or coulomb_c during inversion
     ! Note: inversion%inversion_babc_timescale is later rescaled to SI units (s).
     ! If babc_thck_scale > 0.0, then there is inversion based on a thickness target.
     real(dp) ::  &
          inversion_babc_timescale  = 200.d0,            & !> inversion timescale (yr); must be > 0
          inversion_babc_thck_scale = 100.d0,            & !> thickness inversion scale (m)
          inversion_babc_length_scale = 0.0d0,           & !> diffusive length scale (m) for inversion
          inversion_babc_relax_factor = 0.05d0             !> controls strength of relaxation to default values (unitless)

     ! fields for basin-scale coulomb_c or powerlaw_c inversion
     real(dp), dimension(:,:), pointer ::  &
          grounded_thck_target => null()         !> Observational target for grounded ice thickness

     ! parameters for local deltaT_ocn inversion
     ! Note: deltaT_ocn is in the ocean_data type
     real(dp) ::  &
          inversion_deltaT_ocn_timescale = 200.0d0,      & !> inversion timescale (yr)
          inversion_deltaT_ocn_thck_scale = 100.0d0,     & !> inversion thickness scale (m)
          inversion_deltaT_ocn_temp_scale = 2.0d0,       & !> inversion temperature scale (degC)
          inversion_deltaT_ocn_length_scale = 0.d0         !> diffusive length scale (m) for inversion

     ! fields and parameters for basin-scale deltaT_ocn inversion

     real(dp), dimension(:,:), pointer ::  &
          floating_thck_target => null()         !> Observational target for floating ice thickness

     real(dp) :: &
          inversion_deltaT_basin_relax = 0.0d0,          & !> relax the basin correction toward this value
          inversion_basin_flotation_threshold = 200.d0,  & !> threshold (m) for counting ice as lightly floating/grounded
          inversion_basin_mass_correction = 0.0d0          !> optional mass correction (Gt) for a selected basin

     integer ::  &
          inversion_basin_number_mass_correction = 0       !> integer ID for the basin receiving the correction

  end type glissade_inversion_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_basal_melt_type

     !> Holds fields and parameters relating to basal melting

     ! bmlt fields for grounded and floating ice

     real(dp), dimension(:,:), pointer :: &
          bmlt => null(),                         & !> basal melt rate (> 0 for melt, < 0 for freeze-on)
                                                    !> bmlt = bmlt_ground + bmlt_float (m/s)
          bmlt_applied => null(),                 & !> basal melt rate applied to ice (m/s)
                                                    !>    = 0 for ice-free cells with bmlt > 0
          bmlt_applied_tavg => null(),            & !> basal melt rate applied to ice (m/s, time average)
          bmlt_ground => null(),                  & !> basal melt rate for grounded ice
          bmlt_float => null(),                   & !> basal melt rate for floating ice (m/s)
          bmlt_float_target => null(),            & !> target basal melt rate for floating ice (m/s)
          bmlt_float_external => null(),          & !> external basal melt rate field
          bmlt_float_anomaly => null()              !> basal melt rate anomaly field

     ! masks
     integer, dimension(:,:), pointer :: &
          thermal_forcing_mask => null()            !> = 1 where ice is present, marine-connected and partly floating;
                                                    !> this is where bmlt_float can be nonzero

     real(dp) :: bmlt_float_factor = 1.0d0          !> adjustment factor for external bmlt_float field

     real(dp) :: bmlt_cavity_h0 = 0.0d0             !> scale for reducing melting in sub-shelf cavities (m)
                                                    !> similar to bmlt_float_h0, which is used specifically for MISMIP+

     ! MISMIP+ parameters for Ice1 experiments (BMLT_FLOAT_METHOD_MISMIP)
     ! Note: Parameters with units yr^{-1} are scaled to s^{-1} when used in the code
     real(dp) :: bmlt_float_omega = 0.2d0           !> time scale for basal melting (yr-1)
                                                    !> default value = 0.2 yr^{-1} for MISMIP+ Ice1r
     real(dp) :: bmlt_float_h0 = 75.d0              !> scale for sub-shelf cavity thickness (m)
                                                    !> default value = 75 m for MISMIP+ Ice1r
     real(dp) :: bmlt_float_z0 = -100.d0            !> scale for ice draft, relative to sea level (m)
                                                    !> default value = -100 m for MISMIP+ Ice1r

     ! MISMIP+ parameters for Ice2 experiments (BMLT_FLOAT_METHOD_CONSTANT)
     real(dp) :: bmlt_float_const = 0.d0            !> constant melt rate (m/yr)
                                                    !> set to 100 m/yr for MISMIP+ Ice2r
     real(dp) :: bmlt_float_xlim = 0.d0             !> melting is allowed only for abs(x1) > bmlt_float_xlim
                                                    !> set to 480 km for MISMIP+ Ice2r
     ! parameters for BMLT_FLOAT_METHOD_DEPTH
     ! The sub-shelf melt rate is piecewise linear, generally with greater melting at depth.
     ! The maximum melting and freezing rates are set independently.
     ! The melting/freezing rates fall linearly from the max values to zero over ranges defined by
     !  zmeltmax, zmelt0 and zfrzmax.
     real(dp) :: bmlt_float_depth_meltmax = 10.d0     !> max melt rate at depth (m/yr)
     real(dp) :: bmlt_float_depth_frzmax = 0.d0       !> max freezing rate near surface (m/yr)
     real(dp) :: bmlt_float_depth_zmeltmax = -500.d0  !> depth (m) below which bmlt_float = meltmax
     real(dp) :: bmlt_float_depth_zmelt0 = -200.d0    !> depth (m) where bmlt_float = 0
     real(dp) :: bmlt_float_depth_zfrzmax = -100.d0   !> depth (m) above which bmlt_float = -frzmax

     ! parameters and fields for warm ocean
     ! Where warm_ocean_mask = 1, the sub-shelf melt rate is meltmin above zmeltmin, then increases linearly
     !  to meltmax at zmeltmax. 
     integer, dimension(:,:), pointer :: warm_ocean_mask => null()  !> mask for applying warm ocean melt rate
     real(dp) :: bmlt_float_depth_meltmin = 0.d0      !> min melt rate in warm ocean (m/yr)
     real(dp) :: bmlt_float_depth_zmeltmin = 0.d0     !> depth (m) above which bmlt_float = meltmin

     ! initMIP-Antarctica parameters
     real(dp) :: bmlt_anomaly_tstart = 0.0d0        !> time to start applying the anomaly (yr)
     real(dp) :: bmlt_anomaly_timescale = 0.0d0     !> number of years over which the bmlt_float anomaly is phased in linearly
                                                    !> If set to zero, then the anomaly is applied immediately.


    !> [basal_melt] options (parameters are above)

    integer :: bmlt_float_method = BMLT_FLOAT_METHOD_NONE
    !> basal melt rate for floating ice:
    !> \begin{description}
    !> \item[none] Basal melt rate = 0 for floating ice
    !> \item[mismip] Depth-dependent basal melt rate for floating ice as specified for MISMIP+
    !> \item[constant] Basal melt rate = constant for floating ice (with option to selectively mask out melting)
    !> \item[depth] Depth-dependent basal melt rate for floating ice
    !> \item[file] External basal melt rate field (from input file or coupler)
    !> \item[thermal_forcing] Basal melt rate for floating ice derived from ocean thermal forcing
    !> \end{description}

    integer :: bmlt_float_thermal_forcing_param = BMLT_FLOAT_TF_PARAM_ISMIP6_NONLOCAL
    !> \begin{description}
    !> \item[quadratic] Quadratic parameterization to compute basal melting from thermal forcing
    !> \item[ismip6_local] ISMIP6 local quadratic parameterization to compute basal melting from thermal forcing
    !> \item[ismip6_nonlocal] ISMIP6 nonlocal quadratic parameterization to compute basal melting from thermal forcing
    !> \item[ismip6_nonlocal_slope] ISMIP6 nonlocal quadratic parameterization with slope dependence
    !> \end{description}

    logical :: bmlt_float_init = .false.
    !> if true, then compute bmlt_float at initialization (so it can be written to initial output files)

    integer :: deltaT_ocn_method = DELTAT_OCN_METHOD_NONE
    !> Flag for ocean temperature corrections
    !> \begin{description}
    !> \item[none] deltaT_ocn = 0
    !> \item[inversion] invert for 2D deltaT_ocn to match thickness target
    !> \item[file] read deltaT_ocn from external file
    !> \item[inversion_basin] invert for basin-scale deltaT_ocn
    !> \item[calibrate_basin] calibrate deltaT_ocean to match a basin-scale melt target
    !> \item[dthck_dt] set deltaT_ocn to match a local dH/dt target
    !> \end{description}

    logical :: enable_bmlt_anomaly = .false.
    !> if true, then apply a prescribed anomaly to bmlt_float

    logical :: basal_mbal_in_continuity = .true.
    !> if true, then the basal mass balance is included in the continuity equation

     real(dp) :: bmlt_float_gamma0 = 0.d0                      !> coefficient relating sub-shelf melt rates to thermal forcing (m/yr)

  end type glissade_basal_melt_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_ocean_data_type

     !----------------------------------
     ! fields and parameters for ocean data input
     ! Note: Currently, this type consists of fields for ISMIP6-style ocean input.
     !----------------------------------

     ! ocean grid and basin number
     integer  :: nbasin = 1                         !> number of basins (= 16 for IMBIE2)
     integer  :: nzocn = 1                          !> number of ocean levels
     real(dp) :: dzocn = 0.d0                       !> thickness of ocean levels; nonzero value set in config file

     real(dp), dimension(:), pointer :: &
          zocn => null()                            !> ocean levels (m) where forcing is provided, negative below sea level

     integer, dimension(:), pointer :: &
          basin => null()                           !> basin dimension variable, used for I/O

     real(dp) :: thermal_forcing_basin_min = 0.0d0  !> min value of thermal_forcing_basin (deg K) for nonlocal and nonlocal-slope schemes
     real(dp) :: thermal_forcing_basin_max = 0.0d0  !> max value of thermal_forcing_basin (deg K) for nonlocal and nonlocal-slope schemes
                                                    !> for the default value of 0.0, no upper limit is applied

     ! fields read from input or forcing files

     real(dp), dimension(:,:,:), pointer :: &
          thermal_forcing => null()                 !> 3D thermal forcing (deg K) input to CISM

     real(dp), dimension(:,:), pointer :: &
          thermal_forcing_lsrf => null()            !> 2D thermal forcing (deg K) applied at lower ice surface

     !Note: ocean_data%bmb_float has the same units as climate%smb: kg/m2/yr = mm/yr w.e.
     !      Defined to be positive for melting, since the ISMIP input datasets follow this convention
     real(dp), dimension(:,:), pointer :: &
          bmb_float => null()                       !> 2D basal mass balance (kg/m2/yr) for floating ice

     integer, dimension(:,:), pointer :: &
          basin_number => null()                    !> basin number for each grid cell

     real(dp), dimension(:,:), pointer :: &
          deltaT_ocn => null()                      !> deltaT_ocn in each grid cell (deg K)

     real(dp), dimension(:), pointer :: &
          deltaT_ocn_basin => null()                !> deltaT_ocn applied to each basin (deg K)

     real(dp) :: &
          thermal_forcing_anomaly = 0.0d0,  &       !> thermal forcing anomaly (deg K), applied everywhere
          thermal_forcing_anomaly_tstart = 0.0d0, & !> starting time (yr) for applying or phasing in the anomaly
          thermal_forcing_anomaly_timescale = 0.0d0 !> number of years over which the anomaly is phased in linearly;
                                                    !> if timescale = 0, the full anomaly is applied immediately
     integer :: &
          thermal_forcing_anomaly_basin = 0         !> basin where anomaly is applied;
                                                    !> for default value of 0, apply to all basins


    !> [ocean_data] options (parameters are above)

    integer :: ocean_data_domain = OCEAN_DATA_DOMAIN_FILE
    !> \begin{description}
    !> \item[internal] ocean data computed internally by CISM
    !> \item[file] ocean data read from external file
    !> \item[coupler] ocean data received from coupler via Glad
    !> \end{description}

    logical :: extrapolate_ocean_data = .false.
    !> if true, then CISM extrapolates off-shelf ocean data into shelf cavities

  end type glissade_ocean_data_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_glaciers_type

     !----------------------------------------------------------------
     ! options, fields and parameters for tracking and tuning glaciers
     !----------------------------------------------------------------

     integer :: nglacier = 1                  !> number of glaciers in the global domain

     integer :: ngdiag = 0                    !> CISM index of diagnostic glacier
                                              !> (associated with global cell idiag, jdiag)

     ! inversion options

     integer :: mu_star_method = MU_STAR_METHOD_CONSTANT
     !> \begin{description}
     !> \item[constant] apply spatially uniform mu_star
     !> \item[inversion] invert for glacier-specific mu_star
     !> \item[file] read glacier-specific mu_star from external file
     !> \end{description}

     integer :: alpha_snow_method = ALPHA_SNOW_METHOD_CONSTANT
     !> \begin{description}
     !> \item[constant] apply spatially uniform alpha_snow
     !> \item[inversion] invert for glacier-specific alpha_snow
     !> \item[file] read glacier-specific alpha_snow from external file
     !> \end{description}

     integer :: glacier_powerlaw_c_method = GLACIER_POWERLAW_C_METHOD_CONSTANT
     !> \begin{description}
     !> \item[constant] apply spatially uniform powerlaw_c
     !> \item[inversion] invert for 2D powerlaw_c
     !> \item[file] read 2D powerlaw_c from external file
     !> \end{description}

     ! other options
     integer :: snow_method = SNOW_METHOD_PRECIP_ARTM
     !> \begin{description}
     !> \item[snow] read the snowfall rate directly
     !> \item[precip_artm] compute the snowfall rate from precip and downscaled artm
     !> \end{description}

     logical :: scale_area = .false.
     !> if true, than scale glacier area based on latitude

     logical :: redistribute_advanced_ice = .false.
     !> if true, then thin and redistribute advanced ice in the accumulation zone

     ! parameters
     ! Note: glacier%minthck is currently set at initialization based on model%numerics%thcklim.
     !       glacier%glacier_diag_thck_min is used only for diagnostic area and volume sums;
     !        it does not enter the inversion or dynamics.
     !       Other glacier parameters are declared at the top of module glissade_glacier.
     !       These could be added to the derived type.

     real(dp) :: length_scale_factor = 1.0d0  !> factor used to scale dew and dns;
                                              !> typically equal to the cosine of an average latitude

     real(dp) :: glacier_diag_thck_min = 10.0d0  !> min ice thickness to be included in glacier area and volume diagnostics

     real(dp) :: &
          minthck                             !> min ice thickness (m) to be counted as part of a glacier;
                                              !> currently set based on model%numerics%thcklim

     real(dp) :: &
          glacier_tmlt = -1.d0                        !> spatially uniform temperature threshold for melting (deg C)

     real(dp) :: &
          mu_star_const = 1000.d0,          & ! uniform initial value for mu_star (mm/yr w.e/deg C)
          mu_star_min = 300.d0,             & ! min value of mu_star (mm/yr w.e/deg C)
          mu_star_max = 3000.d0               ! max value of mu_star (mm/yr w.e/deg C)

     real(dp) :: &
          alpha_snow_const = 1.d0,          & ! uniform initial value of alpha_snow (unitless)
          alpha_snow_min = 0.3d0,           & ! min value of alpha_snow
          alpha_snow_max = 3.0d0              ! max value of alpha_snow

     real(dp) ::  &
          !TODO - Change beta_artm_max to 5.0 as in Alps runs?
          beta_artm_max = 3.0,              & ! max magnitude of beta_artm (deg C)
          beta_artm_increment = 0.05d0        ! fixed increment in beta_artm (deg C)

     ! Note: These thresholds assume that artm is a monthly mean, not an instantaneous value
     !       Huss and Hock (2015) have thresholds of 0.5 and 2.5 C
     real(dp) :: &
          glacier_snow_threshold_min = 0.0d0,       & !> air temperature (deg C) below which all precip falls as snow
          glacier_snow_threshold_max = 2.0d0          !> air temperature (deg C) above which all precip falls as rain

     real(dp) :: &
          !TODO - Change baseline date to 1984 as in Alps runs?
          baseline_date = 1980.d0,          & !> baseline date, when glaciers are assumed to be in balance
          rgi_date = 2003.d0,               & !> date of RGI observations
          recent_date = 2010.d0               !> recent date associated with SMB observations for glaciers out of balance

     real(dp) :: &
          thinning_rate_advanced_ice = 0.0d0  !> thinning rate (m/yr) for advanced ice in the accumulation zone;
                                              !> applies when redistribute_advanced_ice = .true.
                                              !> thinned ice volume is redistributed conservatively over the glacier

     real(dp) :: &
          smb_weight_advanced_ice = 1.0d0     !> weight applied to advanced ice in ablation zone during inversion;
                                              !> applied to initially glacier-free cells adjacent to glacier cells
                                              !> typically O(1), with larger values on finer grids
     ! diagnostic scalars

     real(dp) :: &
          total_area = 0.d0,                & !> total area (m^2), summed over all glaciers
          total_volume = 0.d0                 !> total volume (m^3), summed over all glaciers

     integer :: &
          nglacier_active = 0                 !> number of dynamically active glaciers (nonzero area)

     ! 1D arrays with size nglacier

     integer, dimension(:), pointer :: &
          glacierid => null()                 !> glacier ID dimension variable, used for I/O

     ! The following will be allocated with size nglacier, once nglacier is known.

     integer, dimension(:), pointer :: &
          cism_to_rgi_glacier_id => null()    !> maps CISM glacier IDs (1:nglacier) to input RGI glacier IDs

     real(dp), dimension(:), pointer :: &
          area => null(),                   & !> glacier area (m^2)
          volume => null(),                 & !> glacier volume (m^3)
          area_init => null(),              & !> initial glacier area (m^2) based on observations
          volume_init => null(),            & !> initial glacier volume (m^3) based on observations
          area_init_extent => null(),       & !> glacier area (m^2) over the initial ice extent;
                                              !> excludes area where the glacier has advanced
          volume_init_extent => null(),     & !> glacier volume (m^3) over the initial ice extent;
                                              !> excludes volume where the glacier has advanced
          area_target => null(),            & !> glacier area target (m^2) for inversion
          volume_target => null(),          & !> glacier volume target (m^3) for inversion
          mu_star => null(),                & !> glacier-specific parameter relating SMB to monthly mean artm (mm/yr w.e./deg),
                                              !> defined as positive for ablation
          alpha_snow => null(),             & !> glacier-specific multiplicative snow factor (unitless)
          beta_artm => null(),              & !> bias correction to artm (deg C); beta > 0 => increase artm
          smb => null(),                    & !> modeled glacier-average mass balance (mm/yr w.e.)
          smb_obs => null()                   !> observed glacier-average mass balance (mm/yr w.e.), e.g. from Hugonnet et al. (2021)

     ! 2D arrays

     integer, dimension(:,:), pointer :: &
          rgi_glacier_id => null(),         & !> unique glacier ID  based on the Randolph Glacier Inventory
                                              !> first 2 digits give the RGI region;
                                              !> the rest give the number within the region
          cism_glacier_id => null(),        & !> CISM-specific glacier ID, numbered from 1 to nglacier
          cism_glacier_id_init => null(),   & !> cism_glacier_id at initialization, based on rgi_glacier_id
          smb_glacier_id => null(),         & !> integer glacier ID for applying SMB at runtime
          smb_glacier_id_init => null()       !> integer glacier ID for applying SMB;
                                              !> based on cism_glacier_id_init and used for inversion

     real(dp), dimension(:,:), pointer :: &
          dthck_dt_annmean => null(),       & !> annual mean dthck_dt (m/yr)
          snow_annmean => null(),           & !> annual mean snowfall (mm/yr w.e.)
          Tpos_annmean => null(),           & !> annual mean max(artm - tmlt,0) (deg C)
          snow_rgi_annmean => null(),       & !> annual mean snowfall (mm/yr w.e.), RGI date
          Tpos_rgi_annmean => null(),       & !> annual mean max(artm - tmlt,0) (deg C), RGI date
          snow_recent_annmean => null(),    & !> annual mean snowfall (mm/yr w.e.), recent date
          Tpos_recent_annmean => null(),    & !> annual mean max(artm - tmlt,0) (deg C), recent date
          smb_applied_annmean => null()       !> annual mean applied SMB (mm/yr w.e.), = 0 when cell is ice-free

     real(dp), dimension(:,:), pointer :: &
          usrf_target => null(),            & !> target ice surface elevation (m) for the baseline date
          thck_target => null(),            & !> target ice thickness (m) for the baseline date
                                              !> Note: geometry%usrf_obs gives the target for the RGI date
          smb_rgi => null(),                & !> RGI-date SMB field, used for glacier inversion (mm/yr w.e.)
          delta_usrf_rgi => null(),         & !> change in usrf between baseline and RGI climate
          smb_recent => null(),             & !> recent SMB field, including anomaly forcing (mm/yr w.e.)
          delta_usrf_recent => null()         !> change in usrf between baseline and recent climate

     integer, dimension(:,:), pointer :: &
          boundary_mask => null()             !> mask that marks boundary between two glaciers; located at vertices

     integer, dimension(:,:), pointer :: &
          imask => null()                     !> 2D mask; indicates whether glaciers are present in the input file
                                              !> TODO - Remove this field?  Easily derived from initial thickness > 0.


    !> [glaciers] switch

    logical :: enable_glaciers = .false.
    !> if true, then read glacier info at initialization and (optionally)
    !>  tune glacier parameters during the run

  end type glissade_glaciers_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_basal_hydro_type

     !> Holds variables related to basal hydrology
     !> See glissade_bwater.F90 for usage details

     ! hydrology options

     integer :: flux_routing_method = FLUX_ROUTING_METHOD_D8
     !> Flux routing scheme for basal water:
     !> \begin{description}
     !> \item[d8] D8; send flux to lowest-elevation neighbor
     !> \item[dinf] Dinf; divide flux between two lower-elevation neighbors
     !> \item[fd8] FD8; divide flux among all lower-elevation neighbors
     !> \end{description}

     integer :: cavity_open_slide = CAVITY_OPEN_SLIDE_NONE
     !> Treatment of opening by sliding in cavity-sheet equations
     !> \begin{description}
     !> \item[none] no opening by sliding
     !> \item[fixed_ub] opening by sliding with fixed u_b
     !> \item[dynamic_ub] opening by sliding with dynamic u_b
     !> \end{description}

     integer :: cavity_open_melt = CAVITY_OPEN_MELT_NONE
     !> Treatment of opening by melting in cavity-sheet equations
     !> \begin{description}
     !> \item[none] no opening by melting
     !> \item[bmlt] opening by melting based on bmlt_ground
     !> \item[dissipation] opening by melting based on dissipation in cavities
     !> \item[bmlt_dissipation] opening by melting based on bmlt_ground plus dissipation
     !> \end{description}

     !TODO - Add options for closing by creep?
     !       For now, assume closing is proportional to flwa_basal * bwat * effecpress^n

     ! fields related to basal water
     ! Note: Ideally, bwat should have MKS units (m), but currently is scaled.
     !       bwat_diag has units of m; not scaled
     real(dp),dimension(:,:),  pointer :: bwat => null()        !> Local basal water depth
     real(dp),dimension(:,:),  pointer :: bwatflx => null()     !> Basal water flux through a cell per unit cell area (m/s)
     real(dp),dimension(:,:),  pointer :: bmlt_hydro => null()  !> Meltwater source at bed for hydrology (m/s);
                                                                !> sum of local bmlt and englacial/surface meltwater source
     real(dp),dimension(:,:),  pointer :: bwat_diag => null()   !> Steady-state basal water depth diagnosed from bwatflx (m)
     real(dp),dimension(:,:),  pointer :: head => null()        !> Hydraulic head (m)
     real(dp),dimension(:,:),  pointer :: grad_head => null()   !> Gradient of hydraulic head (m/m)

     ! parameters for constant basal water
     ! Note: This parameter applies to the case BWAT_METHOD_CONSTANT.
     real(dp) :: bwat_const = 10.d0              !> constant basal water depth (m)

     ! parameters for local till model
     ! These parameters apply to the case BWAT_METHOD_LOCAL_TILL.
     ! The default values are from Aschwanden et al. (2016) and Bueler and van Pelt (2015).
     real(dp) :: bwat_till_max = 2.0d0           !> maximum water depth in till (m)
     real(dp) :: till_drainage_rate = 1.0d-3     !> uniform drainage rate (m/yr)
     real(dp) :: N_0 = 1000.d0                   !> reference effective pressure (Pa)
     real(dp) :: e_0 = 0.69d0                    !> reference void ratio (dimensionless)
     real(dp) :: C_c = 0.12d0                    !> till compressibility (dimensionless)
                                                 !> Note: The ratio (e_0/C_c) is the key parameter

     ! parameters for steady-state flux-routing model
     !TODO - Add visc_water and omega_hydro? Currently set in glissade_basal_hydro module
     real(dp) :: bwat_source_const = 0.0d0       !> constant melt source at the bed (m/yr)
                                                 !> could be used to represent an englacial or surface source
     real(dp) :: btemp_flow_scale = 0.0d0        !> temperature scale (degC) for transition between thawed and frozen bed;
                                                 !> used to route flow away from cells with a frozen bed;
                                                 !> btemp_flow_scale = 0 => temperature-independent flow
     real(dp) :: btemp_freeze_scale = 0.0d0      !> temperature scale (degC) for transition between thawed and frozen bed;
                                                 !> used to refreeze water beneath cells with a frozen bed;
                                                 !> btemp_freeze_scale = 0 => no refreezing
     ! parameters for macroporous sheet
     real(dp) :: bwat_threshold = 1.0d-3         !> scale over which N ramps down from overburden to a small value (m)

     ! parameters for cavity-sheet
     real(dp) :: bump_height = 0.1d0             !> height scale for bedrock bumps (m);
                                                 !> de Fleurian et al. (2018): 0.1 m (Table 3)
     real(dp) :: bump_wavelength = 2.0d0         !> length scale for bedrock bumps (m); de Fleurian et al. (2018) used 2.0 m
                                                 !> de Fleurian et al. (2018): 2.0 m (Table 3)
     real(dp) :: sliding_speed_fixed = 31.5d0    !> fixed sliding speed (m/yr) for cavity opening by sliding;
                                                 !> de Fleurian et al. (2018): 1.e-6 m/s = 31.5 m/yr (Table 3)
     real(dp) :: flwa_basal = 1.064e-16          !> flow factor A at the bed for cavity closing by creep (Pa^-3 yr-1);
                                                 !> de Fleurian et al. (2018): 3.375e-24 Pa^-3 s^-1 = 1.064e-16 Pa^-3 yr^-1 (Table 3)
     real(dp) :: cavity_closing_const = 0.074    !> constant in the closing equation, v_close = cavity_closing_const * A * bwat * N^n
                                                 !> de Fleurian et al. (2018): 2/(n^n) = 0.074 with n = 3

     ! other parameters for effective pressure
     real(dp) :: effecpress_delta = 0.02d0       !> min value for effecpress N relative to overburden (unitless)
     real(dp) :: bpmp_threshold = 0.1d0          !> temperature range over which N ramps up from a small value to overburden (deg C)
     real(dp) :: haf_threshold = 500.d0          !> thickness at which height above flotation is capped (m)


    !> [basal_hydro] switch

    integer :: bwat_method = BWAT_METHOD_NONE
    !> Basal water depth:
    !> \begin{description}
    !> \item[none] Set to zero everywhere
    !> \item[constant] Set to constant everywhere, to force T = Tpmp
    !> \item[local_till] Local basal till model with constant drainage
    !> \item[flux_routing] Steady-state water routing with flux calculation
    !> \end{description}

  end type glissade_basal_hydro_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_basal_physics_type

     !> Holds variables related to basal physics associated with friction laws.
     !> See glissade_basal_physics.F90 for usage details

     !Note: By default, beta_grounded_min is set to a small nonzero value.
     !      Larger values (~10 to 100 Pa yr/m) might be needed for stability in realistic simulations.
     real(dp) :: beta_grounded_min = 1.0d0  !> minimum value of beta for grounded ice, Pa yr/m (glissade only; scaled during init)
     
     real(dp) :: beta_const = 1000.d0    !> spatially uniform beta for HO dycores, Pa yr/m (scaled during init)
     real(dp) :: beta_small = 1000.d0    !> small beta for sliding over a thawed bed, Pa yr/m (scaled during init)
     real(dp) :: beta_large = 1.0d10     !> large beta to enforce (virtually) no slip, Pa yr/m (scaled during init)

     integer,  dimension(:,:), pointer :: bpmp_mask => null()   !> basal pressure melting point mask; = 1 where Tbed = bpmp, elsewhere = 0
                                                                !> Note: Defined on velocity grid, whereas temp and bpmp are on ice grid

     !WHL TODO - Remove c_space_factor during the config reorg?
     !           Currently used for sliding_law = SLIDING_LAW_COULOMB, *SCHOOF AND *TSAI
     real(dp), dimension(:,:), pointer :: c_space_factor => null()      !> spatial factor for basal shear stress (no dimension)
     real(dp), dimension(:,:), pointer :: c_space_factor_stag => null() !> spatial factor for basal shear stress on staggered grid

     ! parameters for the Zoet-Iverson sliding law
     ! tau_b = N * tan(phi) * [u_b / (u_b + u_t)]^(1/m), Eq. 3 in ZI(2020)
     ! Here, tan(phi) is replaced by coulomb_c
     real(dp) :: zoet_iverson_ut= 200.d0         !> threshold velocity for Zoet-Iverson law (m/yr)
     real(dp) :: zoet_iverson_nmax = 1.0d8       !> max effective pressure for Zoet-Iverson law (Pa)
                                                 !> default value is a high cap, greater than actual overburden pressures

     ! parameters for pseudo-plastic sliding law (based on PISM)
     ! (tau_bx,tau_by) = -tau_c * (u,v) / (u_0^q * |u|^(1-q))
     ! where the yield stress tau_c = tan(phi) * N
     ! N = effective pressure

     real(dp) :: pseudo_plastic_q = 0.5d0        !> exponent for pseudo-plastic law (unitless), 0 <= q <= 1
                                                 !> q = 1 => linear sliding law; q = 0 => plastic; intermediate values => power law
     real(dp) :: pseudo_plastic_u0 = 100.d0      !> threshold velocity for pseudo-plastic law (m/yr)

     ! The following 4 parameters give a linear increase in phi between elevations bedmin and bedmax
     ! Note: These four parameters are used with option SLIDING_LAW_PSEUDO_PLASTIC_OLD
     !       This option will be deprecated but was used for many CESM runs and is kept for backward compatibility
     real(dp) :: pseudo_plastic_phimin =    5.d0 !> min(phi) in pseudo-plastic law, for topg <= bedmin (degrees, 0 < phi < 90)
     real(dp) :: pseudo_plastic_phimax =   40.d0 !> max(phi) in pseudo-plastic law, for topg >= bedmax (degrees, 0 < phi < 90)
     real(dp) :: pseudo_plastic_bedmin = -300.d0 !> bed elevation (m) below which phi = phimin
     real(dp) :: pseudo_plastic_bedmax =  700.d0 !> bed elevation (m) above which phi = phimax

     ! parameters for friction powerlaw
     real(dp) :: friction_powerlaw_k = 8.4d-9    !> coefficient (m y^-1 Pa^-2) for the friction power law based on effective pressure
                                                 !> default value from Bindschadler (1983) based on fits to observations,
                                                 !>  converted to CISM units

     ! Note: powerlaw_c has units of Pa (m/yr)^(-1/powerlaw_m); default value assumes powerlaw_m = 3
     real(dp), dimension(:,:), pointer :: &
          powerlaw_c => null(), &                !> powerlaw_c on staggered grid, Pa (m/yr)^(-1/m)
          coulomb_c => null()                    !> coulomb_c on staggered grid, unitless in range [0,1]

     ! parameters for Weertman-type power law, taub_b = C * u_b^(1/m); used for SLIDING_LAW_SCHOOF AND *_TSAI
     ! The default values are from Asay-Davis et al. (2016).
     ! The value of powerlaw_c suggested by Tsai et al. (2015) is 7.624d6 Pa m^(-1/3) s^(1/3).
     ! This value can be converted to CISM units by dividing by scyr^(1/3), to obtain 2.413d4 Pa m^(-1/3) yr^(1/3).
     ! Note: The Tsai et al. Coulomb friction law uses coulomb_c above, with
     !       effective pressure N as in Leguy et al. (2014) with p_ocean_penetration = 1.

     real(dp) :: powerlaw_c_const = 1.0d4        !> friction coefficient in power law, units of Pa m^(-1/3) yr^(1/3)
     real(dp) :: powerlaw_m = 3.d0               !> exponent in power law (unitless)
     real(dp) :: powerlaw_c_max = 1.0d5          !> max value of powerlaw_c, Pa (m/yr)^(-1/3)
     real(dp) :: powerlaw_c_min = 1.0d2          !> min value of powerlaw_c, Pa (m/yr)^(-1/3)

     ! parameters for Coulomb basal friction law
     ! Notes: coulomb_c_max = 1.0 to cap effecpress at overburden
     !        The appropriate value of coulomb_c_min can depend on how much N is reduced below overburden.
     ! TODO: Reduce coulomb_c_const to 0.1?
     real(dp) :: coulomb_c_const = 0.42d0        !> basal stress constant; unitless in range [0,1]
     real(dp) :: coulomb_c_max = 1.0d0           !> max value of coulomb_c, unitless
     real(dp) :: coulomb_c_min = 1.0d-3          !> min value of coulomb_c, unitless

     ! parameters for Schoof basal friction law
     ! Typically, schoof_n (from Eq. 2 in Schoof 2005) has the same value as powerlaw_m, but this is not required.
     ! The parameters gamma and p are not in Schoof (2005) but can be used to relate Cc and Cp during inversion,
     !  with either the Schoof law or the Tsai law, if powerlaw_c_method = POWERLAW_C_METHOD_FROM_COULOMB_C.
     !  The default values here are consistent with Cc_max = 1.0, Cp_max = 1.e5, Cc_const = 1.0, Cp_const ~ 2.e4
     real(dp) :: schoof_n = 3.0d0                !> exponent in the Schoof basal friction law;
                                                 !> modulates the transition between powerlaw and coulomb behavior
     real(dp) :: schoof_gamma = 1.0d5            !> parameters in the relation Cp = gamma * Cc^p;
     real(dp) :: schoof_p = 0.7d0                !> only used if nonzero values are set in the config file

     ! The next four parameters apply when elevation_based_coulomb_c = .true.
     real(dp) :: coulomb_c_hi =  0.50d0          !> coulomb_c value at high bed elevation, topg >= bed_hi
     real(dp) :: coulomb_c_lo =  0.05d0          !> coulomb_c value at low bed elevation, topg <= bed_lo
     real(dp) :: coulomb_c_bed_hi =    0.d0      !> bed elevation (m) above which coulomb_c = coulomb_c_const_hi
     real(dp) :: coulomb_c_bed_lo = -500.d0      !> bed elevation (m) below which coulomb_c = coulomb_c_const_lo

     ! parameters for older form of Coulomb friction sliding law (default values from Pimentel et al. 2010)
     ! Pimentel et al. have coulomb_c = 0.84*m_max, where m_max = coulomb_bump_max_slope
     !TODO - Remove these constants? They are used for babc option SLIDING_LAW_COULOMB,
     !       which is functionally equivalent to the Schoof law.
     real(dp) :: coulomb_bump_wavelength = 2.0d0 !> bedrock wavelength at subgrid scale precision (m)
     real(dp) :: flwa_basal_const = 1.0d-16      !> flow factor at the bed for sliding_law = 10 (Pa^{-n} yr^{-1});
                                                 !> 1.0e-16 Pa^-3 yr^-1 = 3.1688e-24 Pa^-3 s^-1, as in MISMIP3d (Leguy et al.)
     real(dp) :: coulomb_bump_max_slope = 0.5d0  !> maximum bed bump slope at subgrid scale precision (no dimension)

     ! parameter to limit the min value of beta for various power laws
     real(dp) :: beta_powerlaw_umax = 0.0d0      !> upper limit of ice speed (m/yr) when evaluating powerlaw beta
                                                 !> Where u > umax, let u = umax when evaluating beta(u)

     ! Note: A basal process model is not currently supported, but a specified mintauf can be passed to subroutine calcbeta
     !       to simulate a plastic bed.
     real(dp),dimension(:,:), pointer :: mintauf => null() ! Bed strength (yield stress) calculated with basal process model

     !TODO - At some point, move the remaining fields to the basal_hydro type
     !       Moving them will break some existing config files.
     ! fields related to the effective pressure
     real(dp), dimension(:,:), pointer :: effecpress => null()          !> effective pressure (Pa)
     real(dp), dimension(:,:), pointer :: effecpress_stag => null()     !> effective pressure on staggered grid (Pa)

     ! parameters for reducing the effective pressure where the bed is connected to the ocean
     real(dp) :: p_ocean_penetration = 0.0d0           !> p-exponent for ocean penetration; N weighted by (1-Hf/H)^p (0 <= p <= 1)


    !> [basal_physics] keys (others are above)

    integer :: sliding_law = SLIDING_LAW_NO_SLIP

    logical :: use_c_space_factor = .false.
    !> Flag that indicates whether to multiply beta by a 2D scalar field

    integer :: beta_limit_method = BETA_LIMIT_METHOD_ABSOLUTE
    !> Flag for setting minimum beta beneath grounded ice
    !> \begin{description}
    !> \item[absolute] absolute limit given by beta_grounded_min
    !> \item[fground_scaled] limited using beta_grounded_min, then multiplied by f_ground
    !> \end{description}

    integer :: powerlaw_c_method = POWERLAW_C_METHOD_CONSTANT
    !> Flag for basal powerlaw_c options
    !> \begin{description}
    !> \item[constant] powerlaw_c = spatially uniform constant
    !> \item[inversion] invert for 2D powerlaw_c
    !> \item[file] read 2D powerlaw_c from external file
    !> \item[inversion_basin] invert for basin-scale powerlaw_c
    !> \item[from_coulomb_c] relate to coulomb_c: Cp = gamma*Cc^p
    !> \end{description}

    integer :: coulomb_c_method = COULOMB_C_METHOD_CONSTANT
    !> Flag for basal coulomb_c options
    !> \begin{description}
    !> \item[constant] coulomb_c = spatially uniform constant
    !> \item[inversion] invert for 2D coulomb_c
    !> \item[file] read 2D coulomb_c from external file
    !> \item[inversion_basin] invert for basin-scale coulomb_c
    !> \end{description}

    logical :: elevation_based_coulomb_c = .false.
    !> Flag that indicates whether coulomb_c depends on elevation

    integer :: effecpress_method = EFFECPRESS_METHOD_OVERBURDEN
    !> Flag that describes effective pressure calculation for HO dyn core: 
    !> \begin{description}
    !> \item[overburden] N = overburden pressure, rhoi*grav*thck
    !> \item[bpmp] N is proportional to height above flotation
    !> \item[bwat] N is reduced based on basal water thickness from active hydrology
    !> \item[cavity_sheet] N depends on cavity opening and closing
    !> \item[bwat_bvp] N is reduced based on basal water thickness, following Bueler/van Pelt
    !> \item[5] N is reduced where the bed is at or near the pressure melting point
    !> \end{description}

    integer :: btrc_method = BTRC_METHOD_NONE
    !> Basal slip coefficient (used for glide and for glissade local SIA; not glissade HO):
    !> \begin{description}
    !> \item[none] Set equal to zero everywhere
    !> \item[constant] Set to (non--zero) constant
    !> \item[constant_bwat] Set to (non--zero) constant where basal water is present, otherwise to zero
    !> \item[constant_bpmp] Set to (non--zero) constant where bed temperature is at pressure melting point, otherwise to zero
    !> \item[4] linear function of basal melt rate
    !> \item[5] \texttt{tanh} function of basal water depth 
    !> \end{description}

    real(dp) :: btrc_const = 0.d0      ! m yr^{-1} Pa^{-1} (gets scaled during init)

  end type glissade_basal_physics_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_isos_elastic_type
     !> Holds data used by isostatic adjustment calculations

     real(dp) :: d = 0.24d25                !> flexural rigidity (N m)
     real(dp) :: lr                         !> radius of relative stiffness
     real(dp) :: a                          !> radius of disk
     real(dp) :: c1,c2,cd3,cd4              !> coefficients
     real(dp), dimension(:,:), pointer :: w => null() !> matrix operator for lithosphere deformation
     integer :: wsize                       !> size of operator (0:rbel_wsize, 0:rbel_wsize), operator is axis symmetric
  end type glissade_isos_elastic_type

  type glissade_isostasy_type
     !> contains isostasy configuration

     integer :: lithosphere_method = LITHOSPHERE_METHOD_ELASTIC
     !> method for calculating equilibrium bedrock depression
     !> \begin{description}
     !> \item[local] local lithosphere, equilibrium bedrock depression is found using Archimedes' principle
     !> \item[elastic] elastic lithosphere, flexural rigidity is taken into account
     !> \end{description}

     integer :: asthenosphere_method = ASTHENOSPHERE_METHOD_RELAXING
     !> method for approximating the mantle
     !> \begin{description}
     !> \item[fluid] fluid mantle, isostatic adjustment happens instantaneously
     !> \item[relaxing] relaxing mantle, exp adjust toward (relx - load), constant timescale
     !> \item[relaxing_latvar] relaxing mantle, exp adjust toward (relx - load), laterally varying timescale
     !> \end{description}

    integer :: input_topg_method = INPUT_TOPG_METHOD_SEPARATE_RELX

    !> relaxed topography:
    !> \begin{description}
    !> \item[separate_relx] get current topo (topg) and relaxed topo (relx) from separate input fields
    !> \item[relx_from_topg] first time slice of input topo is relaxed
    !> \item[relx_from_load] first time slice of input topo is in isostatic equilibrium
    !> \end{description}

     real(dp) :: tau_relax_const = 3000.d0    ! characteristic time constant of relaxing mantle (yr)
     real(dp) :: load_update_interval = 10.d0 ! update interval for elastic load (yr)
     integer :: nlith                         ! update the load every nlith time steps; tinc * nlith = load_update_interval
     logical :: new_load = .false.            ! set to true if there is a new surface load
     type(glissade_isos_elastic_type) :: rbel               ! structure holding elastic lithosphere setup

     real(dp),dimension(:,:),pointer :: relx => null()  ! elevation of relaxed topography (m)
     real(dp),dimension(:,:),pointer :: load => null()  ! deflection due to applied load on lithosphere (m)
                                                        ! defined as positive for downward deflection
     real(dp),dimension(:,:),pointer :: load_factors => null() ! temporary used for load calculation
     real(dp),dimension(:,:),pointer :: tau_relax => null()    ! 2D mantle relaxation timescale (yr)


    !> [isostasy] switch and lithosphere parameter

    logical :: compute_isostasy = .false.
    !> if true, then compute isostatic adjustment of the bed

    real(dp) :: flexural_rigidity = 0.24d25  !> flexural rigidity (N m) of the elastic lithosphere; copied to rbel%d

  end type glissade_isostasy_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_time_type

    !> [time] keys and the time-stepping state
    real(dp) :: tstart =     0.d0 !> starting time
                                  !> Note: For restarts, tstart is set to the time the restart begins
    real(dp) :: tend   =  1000.d0 !> end time
    real(dp) :: time   =     0.d0 !> main loop counter in years
    real(dp) :: tinc   =     1.d0 !> time step of main loop in years (key dt)
    real(dp) :: dt     =     0.d0     !> ice dynamics timestep (s)
    integer  :: nsteps_per_year = 0   !> number of time steps per year (for dt_input_option = nsteps_per_year)
    integer  :: tstep_count = 0   !> number of time steps since the start of the simulation

    integer :: dt_input_option = DT_INPUT_OPTION_DT
    !> \begin{description}
    !> \item[dt] Input dt in years
    !> \item[nsteps_per_year] Input dt in steps per year
    !> \end{description}

    ! CFL diagnostics and adaptive subcycling
    real(dp) :: adv_cfl_dt = 0.0d0  ! maximum allowable dt (yrs) based on advective CFL (calculated by model for each time step)
    real(dp) :: diff_cfl_dt = 0.0d0 ! maximum allowable dt (yrs) based on diffusive CFL (calculated by model for each time step)
    real(dp) :: adaptive_cfl_threshold = 0.0d0  ! adaptively subcycle the advection when advective CFL exceeds this value
                                                ! (zero value means no adaptive subcycling)

  end type glissade_time_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_diagnostics_type

    !> [diagnostics] keys and the local diagnostic point
    real(dp) :: dt_diag = 0.d0     ! diagnostic interval (write diagnostic output every dt_diag years)
                                   ! dt_diag = 0 => never write diagnostic output
    integer  :: ndiag = 0          ! diagnostic interval (write diagnostic output every ndiag timesteps)
                                   ! ndiag = 0 => never write diagnostic output
    integer  :: idiag = 1          ! global grid indices for diagnostic point
    integer  :: jdiag = 1          ! 
    integer  :: idiag_local = 1    ! local grid indices for diagnostic point
    integer  :: jdiag_local = 1
    integer  :: rdiag_local = 0    ! task number for diagnostic point

    integer :: dm_dt_diag = DM_DT_DIAG_KG_PER_S
    !> \begin{description}
    !> \item[kg_per_s] Write dmass/dt diagnostic in units of kg/s
    !> \item[gt_per_yr] Write dmass/dt diagnostic in units of Gt/yr
    !> \end{description}

    integer :: diag_ice_threshold = DIAG_ICE_THRESHOLD_THCKLIM
    !> \begin{description}
    !> \item[zero] Include cells with H > 0 in global diagnostics
    !> \item[thcklim] Include cells with H > thcklim in global diagnostics
    !> \end{description}

  end type glissade_diagnostics_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_initialization_type

    !> [initialization] keys
    logical :: adjust_input_thickness = .false.
    !> if true, then adjust thck to maintain usrf, instead of deriving usrf from topg and thck

    logical :: smooth_input_usrf = .false.
    !> if true, then apply Laplacian smoothing to usrf at initialization

  end type glissade_initialization_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_numerics_type

    !> [numerics] keys
    real(dp) :: thcklim =   100.d0 !> min thickness for computing ice dynamics (m) 
    real(dp) :: thcklim_temp = 1.d0    !> min thickness for computing vertical temperature (m) (higher-order only)
    real(dp) :: thck_gradient_ramp = 0.d0 !> thickness scale over which gradients increase from zero to full value (HO only)
    real(dp) :: periodic_offset_ew = 0.d0 !> optional periodic_offsets for ismip-hom and similar tests;
    real(dp) :: periodic_offset_ns = 0.d0 !> may be needed to ensure continuous ice geometry at edges of the global domain

    integer :: compute_blocks = COMPUTE_BLOCKS_ALL
    !> \begin{description}
    !> \item[all] Compute on all blocks in the global domain; one task per block
    !> \item[active_only] Compute on active blocks only; one task per active block, no task for some or all inactive blocks
    !> \item[inquire] Inquire how many active blocks there are, in prep for resubmitting with active_only
    !> \end{description}

    logical :: reproducible_sums = .false.
    !> if true, then compute reproducible global sums
    !> (independent of the number of tasks)

  end type glissade_numerics_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glissade_global_type    ! type containing all of the above for an ice sheet model instance
    integer              :: model_id = 0 !> identifier of this model instance (1..N, where N is the number of ice sheets in this run)
    type(glissade_grid_type)     :: grid
    type(glissade_options_type)  :: options
    type(glissade_geometry_type) :: geometry
    type(glissade_geomderv_type) :: geomderv
    type(glissade_velocity_type) :: velocity
    type(glissade_velocity_solver_type) :: velocity_solver
    type(glissade_advection_type) :: advection
    type(glissade_scalars_type)  :: scalars
    type(glissade_mass_flux_type):: mass_flux
    type(glissade_stress_type) :: stress   
    type(glissade_climate_type)  :: climate
    type(glissade_calving_type)  :: calving
    type(glissade_removal_type)  :: removal
    type(glissade_lateral_melt_type) :: lateral_melt
    type(glissade_temperature_type)   :: temperature
    type(glissade_geothermal_type) :: geothermal
    type(glissade_flow_law_type) :: flow_law
    type(glissade_basal_hydro_type)  :: basal_hydro
    type(glissade_basal_physics_type):: basal_physics
    type(glissade_grounding_line_type) :: grounding_line
    type(glissade_basal_melt_type)   :: basal_melt
    type(glissade_ocean_data_type)   :: ocean_data
    type(glissade_glaciers_type)  :: glaciers
    type(glissade_inversion_type):: inversion
    type(glimmer_funits) :: funits
    type(glissade_time_type)     :: time
    type(glissade_diagnostics_type) :: diagnostics
    type(glissade_initialization_type) :: initialization
    type(glissade_numerics_type) :: numerics
    type(glimmap_proj)   :: projection
    type(profile_type)   :: profile
    type(glissade_isostasy_type)  :: isostasy
    type(parallel_type)  :: parallel

  end type glissade_global_type

contains

  subroutine glissade_allocarr(model)    

    !> Allocates the model arrays, and initialises some of them to zero.

    use glimmer_log
    use glimmer_coordinates, only: coordsystem_allocate

    implicit none

    type(glissade_global_type),intent(inout) :: model

    integer :: ewn,nsn,upn               !> local array dimensions
    integer :: global_ewn, global_nsn    !> global array dimensions
    integer :: n

    ! for simplicity, copy these values...
    ewn = model%grid%ewn
    nsn = model%grid%nsn
    upn = model%grid%upn
    global_ewn = model%parallel%global_ewn
    global_nsn = model%parallel%global_nsn
    
    ! horizontal coordinates

    allocate(model%grid%x0(ewn-1)); model%grid%x0 = 0.d0  ! velocity grid
    allocate(model%grid%y0(nsn-1)); model%grid%y0 = 0.d0
    allocate(model%grid%x1(ewn)); model%grid%x1 = 0.d0    ! ice grid (for scalars)
    allocate(model%grid%y1(nsn)); model%grid%y1 = 0.d0

    allocate(model%grid%x0_global(global_ewn-1)); model%grid%x0_global = 0.d0  ! velocity grid
    allocate(model%grid%y0_global(global_nsn-1)); model%grid%y0_global = 0.d0
    allocate(model%grid%x1_global(global_ewn)); model%grid%x1_global = 0.d0    ! ice grid (for scalars)
    allocate(model%grid%y1_global(global_nsn)); model%grid%y1_global = 0.d0

    ! vertical sigma coordinates
    ! If we already have sigma, don't reallocate

    if (associated(model%grid%sigma)) then
       if (size(model%grid%sigma) /= upn) then
          call write_log('Wrong number of sigma levels given',GM_FATAL)
       end if
    else
       allocate(model%grid%sigma(upn))
    endif

    allocate(model%grid%stagsigma(upn-1))
    allocate(model%grid%stagwbndsigma(0:upn))  !MJH added (0:upn) as separate variable

    ! latitude and longitude
    call coordsystem_allocate(model%grid%ice_grid, model%grid%lat)
    call coordsystem_allocate(model%grid%ice_grid, model%grid%lon)

    ! ice domain mask (to identify active blocks)
    call coordsystem_allocate(model%grid%ice_grid, model%grid%ice_domain_mask)

    ! mask to identify cells at global domain edge
    call coordsystem_allocate(model%grid%ice_grid, model%grid%global_edge_mask)

    ! temperature arrays

    !NOTE: The temperature and flow factor live on
    !       the staggered vertical grid, with temp and flwa defined at the
    !       center of each layer k = 1:upn-1.  The temperature (but not flwa)
    !       is defined at the upper surface (k = 0) and lower surface (k = upn).

    allocate(model%temperature%temp(0:upn,1:ewn,1:nsn))
    ! tempunstag has the same horizontal grid as temp, but is defined
    ! on the unstaggered vertical grid (levels 1:upn)
    allocate(model%temperature%tempunstag(upn,1:ewn,1:nsn))
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%temperature%flwa)
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%temperature%dissip)

    ! MJH - Set temp and flwa to physically unrealistic values so we can tell later if 
    !       arrays were initialized correctly
    model%temperature%temp(:,:,:) = unphys_val  ! large negative number
    model%temperature%flwa(:,:,:) = unphys_val
    model%temperature%dissip(:,:,:) = 0.d0
    if (associated(model%temperature%tempunstag)) &
         model%temperature%tempunstag(:,:,:) = unphys_val

    call coordsystem_allocate(model%grid%ice_grid,  model%temperature%bheatflx)
    call coordsystem_allocate(model%grid%ice_grid,  model%temperature%bpmp)
    call coordsystem_allocate(model%grid%ice_grid,  model%temperature%btemp)
    call coordsystem_allocate(model%grid%ice_grid,  model%temperature%btemp_ground)
    call coordsystem_allocate(model%grid%ice_grid,  model%temperature%btemp_float)
    call coordsystem_allocate(model%grid%ice_grid,  model%temperature%ucondflx)
    call coordsystem_allocate(model%grid%ice_grid,  model%temperature%flow_enhancement_factor)

    call coordsystem_allocate(model%grid%ice_grid,  model%basal_hydro%bwat)
    call coordsystem_allocate(model%grid%ice_grid,  model%basal_hydro%bwatflx)
    if (model%basal_hydro%bwat_method == BWAT_METHOD_FLUX_ROUTING) then
       call coordsystem_allocate(model%grid%ice_grid,  model%basal_hydro%bmlt_hydro)
       call coordsystem_allocate(model%grid%ice_grid,  model%basal_hydro%bwat_diag)
       call coordsystem_allocate(model%grid%ice_grid,  model%basal_hydro%head)
       call coordsystem_allocate(model%grid%ice_grid,  model%basal_hydro%grad_head)
    endif

    call coordsystem_allocate(model%grid%ice_grid, model%temperature%bfricflx)
    call coordsystem_allocate(model%grid%ice_grid, model%temperature%bhydroflx)
    ! water fraction and enthalpy live at the midpoint of each layer (with temp and flwa)
    ! enthalpy (like temp) is defined at the upper and lower surfaces as well
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%temperature%waterfrac)
    allocate(model%temperature%enthalpy(0:upn,1:ewn,1:nsn))
    model%temperature%enthalpy(:,:,:) = 0.d0

    ! velocity arrays

    call coordsystem_allocate(model%grid%velo_grid, upn, model%velocity%uvel)
    call coordsystem_allocate(model%grid%velo_grid, upn, model%velocity%vvel)
    call coordsystem_allocate(model%grid%velo_grid, upn, model%velocity%velnorm)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%btrc)
    call coordsystem_allocate(model%grid%velo_grid, upn, model%velocity%resid_u)
    call coordsystem_allocate(model%grid%velo_grid, upn, model%velocity%resid_v)
    call coordsystem_allocate(model%grid%velo_grid, upn, model%velocity%rhs_u)
    call coordsystem_allocate(model%grid%velo_grid, upn, model%velocity%rhs_v)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%uvel_2d)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%vvel_2d)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%uvel_mean)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%vvel_mean)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%velnorm_mean)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%usfc_obs)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%vsfc_obs)
    call coordsystem_allocate(model%grid%ice_grid,  upn, model%velocity%wvel)

    ! The following are on the extended staggered grid, which is the same size as the ice grid.
    call coordsystem_allocate(model%grid%ice_grid,  upn, model%velocity%uvel_extend)
    call coordsystem_allocate(model%grid%ice_grid,  upn, model%velocity%vvel_extend)
    call coordsystem_allocate(model%grid%ice_grid,  model%velocity%uvel_2d_extend)
    call coordsystem_allocate(model%grid%ice_grid,  model%velocity%vvel_2d_extend)

    call coordsystem_allocate(model%grid%velo_grid, model%velocity%beta)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%beta_internal)
    call coordsystem_allocate(model%grid%ice_grid, model%velocity%unstagbeta)
    ! Set beta and unstagbeta to physically unrealistic values so we can tell later 
    ! if these fields were read correctly from an input file
    model%velocity%beta(:,:) = unphys_val   ! large negative number
    model%velocity%unstagbeta(:,:) = unphys_val

    call coordsystem_allocate(model%grid%velo_grid, model%velocity%kinbcmask)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%umask_no_penetration)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%vmask_no_penetration)
    call coordsystem_allocate(model%grid%ice_grid, model%velocity%divu)
    call coordsystem_allocate(model%grid%ice_grid, model%velocity%shear)

    ! higher-order stress arrays

    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%stress%efvs)
    call coordsystem_allocate(model%grid%ice_grid, model%stress%efvs_vertavg)
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%stress%tau%scalar) 
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%stress%tau%xz)
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%stress%tau%yz)
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%stress%tau%xx)
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%stress%tau%yy)
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%stress%tau%xy)
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%velocity%strain_rate%scalar) 
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%velocity%strain_rate%xz)
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%velocity%strain_rate%yz)
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%velocity%strain_rate%xx)
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%velocity%strain_rate%yy)
    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%velocity%strain_rate%xy)
    call coordsystem_allocate(model%grid%velo_grid, model%stress%btractx)
    call coordsystem_allocate(model%grid%velo_grid, model%stress%btracty)
    call coordsystem_allocate(model%grid%velo_grid, model%stress%btract)
    call coordsystem_allocate(model%grid%ice_grid, model%stress%btractx_extend)
    call coordsystem_allocate(model%grid%ice_grid, model%stress%btracty_extend)
    call coordsystem_allocate(model%grid%velo_grid, model%stress%taudx)
    call coordsystem_allocate(model%grid%velo_grid, model%stress%taudy)

    ! geometry arrays
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%thck)
    call coordsystem_allocate(model%grid%velo_grid, model%geometry%stagthck)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%usrf)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%lsrf)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%topg)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%usrf_obs)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%dthck_dt)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%dthck_dt_obs)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%dthck_dt_obs_basin)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%cell_area)

    !Note: model%geometry%tracers and related arrays are allocated later, in glissade_transport_setup_tracers,
    !      once we know model%geometry%ntracers

    call coordsystem_allocate(model%grid%velo_grid, model%geomderv%dthckdew)
    call coordsystem_allocate(model%grid%velo_grid, model%geomderv%dthckdns)
    call coordsystem_allocate(model%grid%velo_grid, model%geomderv%dusrfdew)
    call coordsystem_allocate(model%grid%velo_grid, model%geomderv%dusrfdns)

    call coordsystem_allocate(model%grid%ice_grid, model%geometry%ice_mask)
    call coordsystem_allocate(model%grid%velo_grid, model%geometry%ice_mask_stag)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%floating_mask)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%grounded_mask)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%marine_connection_mask)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%ice_sheet_mask)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%ice_cap_mask)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%ice_fraction_retreat_mask)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%removal_thck)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%removal_rate)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%removal_rate_tavg)
    call coordsystem_allocate(model%grid%ice_grid, model%geometry%reference_thck)

    ! mass flux arrays
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%sfc_mbal_flux)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%sfc_mbal_flux_tavg)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%basal_mbal_flux)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%basal_mbal_flux_tavg)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%calving_flux)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%calving_flux_tavg)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%latmelt_flux)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%latmelt_flux_tavg)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%removal_flux)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%removal_flux_tavg)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%gl_flux_east)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%gl_flux_north)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%gl_flux)
    call coordsystem_allocate(model%grid%ice_grid, model%mass_flux%gl_flux_tavg)

    call coordsystem_allocate(model%grid%ice_grid, upn-1, model%geometry%ice_age)
    call coordsystem_allocate(model%grid%ice_grid,  model%geometry%thck_old)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%velo_sfc_obs)
    call coordsystem_allocate(model%grid%velo_grid, model%velocity%velo_sfc)
    call coordsystem_allocate(model%grid%ice_grid,  model%geometry%f_flotation)
    call coordsystem_allocate(model%grid%velo_grid, model%geometry%f_ground)
    call coordsystem_allocate(model%grid%ice_grid,  model%geometry%f_ground_cell)

    ! Basal Physics
    call coordsystem_allocate(model%grid%velo_grid, model%basal_physics%bpmp_mask)
    call coordsystem_allocate(model%grid%ice_grid, model%basal_physics%effecpress)
    call coordsystem_allocate(model%grid%velo_grid, model%basal_physics%effecpress_stag)
    call coordsystem_allocate(model%grid%ice_grid, model%basal_physics%c_space_factor)
    call coordsystem_allocate(model%grid%velo_grid, model%basal_physics%c_space_factor_stag)
    call coordsystem_allocate(model%grid%velo_grid, model%basal_physics%mintauf)

    ! bmlt arrays
    call coordsystem_allocate(model%grid%ice_grid,  model%basal_melt%bmlt)
    call coordsystem_allocate(model%grid%ice_grid,  model%basal_melt%bmlt_applied)
    call coordsystem_allocate(model%grid%ice_grid,  model%basal_melt%bmlt_applied_tavg)
    call coordsystem_allocate(model%grid%ice_grid,  model%basal_melt%bmlt_ground)

    call coordsystem_allocate(model%grid%ice_grid, model%ocean_data%basin_number)
    call coordsystem_allocate(model%grid%ice_grid, model%basal_melt%bmlt_float)
    call coordsystem_allocate(model%grid%ice_grid, model%basal_melt%bmlt_float_anomaly)
    call coordsystem_allocate(model%grid%ice_grid, model%basal_melt%warm_ocean_mask)
    call coordsystem_allocate(model%grid%ice_grid, model%basal_melt%bmlt_float_target)
    call coordsystem_allocate(model%grid%ice_grid, model%basal_melt%bmlt_float_external)
    call coordsystem_allocate(model%grid%ice_grid, model%basal_melt%thermal_forcing_mask)
    ! Note: nzocn and nbasin should be set in the [ocean_data] section of the config file
    ! The ocean arrays are allocated in every run, since they are used by several options (basal melt from
    !  thermal forcing, deltaT_ocn inversion, lateral melt, GLAD coupling). With the defaults (nzocn = 1,
    !  nbasin = 1), this costs only a few 2D fields in runs that do not use them.
    call coordsystem_allocate(model%grid%ice_grid, model%ocean_data%nzocn, &
                              model%ocean_data%thermal_forcing)
    call coordsystem_allocate(model%grid%ice_grid, model%ocean_data%thermal_forcing_lsrf)
    call coordsystem_allocate(model%grid%ice_grid, model%ocean_data%bmb_float)
    call coordsystem_allocate(model%grid%ice_grid, model%ocean_data%deltaT_ocn)
    allocate(model%ocean_data%deltaT_ocn_basin(model%ocean_data%nbasin))
    model%ocean_data%deltaT_ocn_basin = 0.0d0   ! computed only for some deltaT_ocn options, but always written

    ! glacier options
    if (model%glaciers%enable_glaciers) then
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%rgi_glacier_id)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%cism_glacier_id)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%cism_glacier_id_init)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%smb_glacier_id)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%smb_glacier_id_init)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%dthck_dt_annmean)
       call coordsystem_allocate(model%grid%velo_grid, model%glaciers%boundary_mask)

       ! Note: The recent and RGI fields are used for glacier inversion
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%usrf_target)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%thck_target)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%smb_rgi)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%delta_usrf_rgi)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%smb_recent)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%delta_usrf_recent)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%snow_annmean)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%Tpos_annmean)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%snow_rgi_annmean)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%Tpos_rgi_annmean)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%snow_recent_annmean)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%Tpos_recent_annmean)
       call coordsystem_allocate(model%grid%ice_grid, model%glaciers%smb_applied_annmean)

       ! Allocate arrays with dimension(nglacier)
       ! Note: nglacier = 1 by default, but can be changed in subroutine glissade_glacier_init
       !        after reading the input file.  If so, these arrays will be reallocated.
       !       On restart, nglacier is read from the restart file before calling glissade_allocarr,
       !        so these allocations will be correct.
       allocate(model%glaciers%glacierid(model%glaciers%nglacier))
       allocate(model%glaciers%cism_to_rgi_glacier_id(model%glaciers%nglacier))
       allocate(model%glaciers%area(model%glaciers%nglacier))
       allocate(model%glaciers%volume(model%glaciers%nglacier))
       allocate(model%glaciers%area_init(model%glaciers%nglacier))
       allocate(model%glaciers%volume_init(model%glaciers%nglacier))
       allocate(model%glaciers%area_init_extent(model%glaciers%nglacier))
       allocate(model%glaciers%volume_init_extent(model%glaciers%nglacier))
       allocate(model%glaciers%area_target(model%glaciers%nglacier))
       allocate(model%glaciers%volume_target(model%glaciers%nglacier))
       allocate(model%glaciers%mu_star(model%glaciers%nglacier))
       allocate(model%glaciers%alpha_snow(model%glaciers%nglacier))
       allocate(model%glaciers%beta_artm(model%glaciers%nglacier))
       allocate(model%glaciers%smb(model%glaciers%nglacier))
       allocate(model%glaciers%smb_obs(model%glaciers%nglacier))
    endif

    ! inversion and basal physics arrays
    call coordsystem_allocate(model%grid%velo_grid,model%basal_physics%powerlaw_c)
    call coordsystem_allocate(model%grid%velo_grid,model%basal_physics%coulomb_c)

    if (model%basal_physics%coulomb_c_method == COULOMB_C_METHOD_INVERSION_BASIN .or. &
        model%basal_physics%powerlaw_c_method == POWERLAW_C_METHOD_INVERSION_BASIN) then
       if (model%ocean_data%nbasin < 1) then
          call write_log ('Must set nbasin >= 1 for basin-scale inversion of C_c or C_p', GM_FATAL)
       endif
       call coordsystem_allocate(model%grid%ice_grid, model%inversion%grounded_thck_target)
    endif

    if (model%basal_melt%deltaT_ocn_method == DELTAT_OCN_METHOD_INVERSION_BASIN) then
       if (model%ocean_data%nbasin < 1) then
          call write_log ('Must set nbasin >= 1 for basin-scale deltaT_ocn inversion', GM_FATAL)
       endif
       call coordsystem_allocate(model%grid%ice_grid, model%inversion%floating_thck_target)
    endif

    ! basin diagnostic arrays
    ! Note: These are allocated even if there is just one basin,
    !       but they are intended for use with multiple basins.
    if (model%ocean_data%nbasin >= 1) then
       allocate(model%ocean_data%basin(model%ocean_data%nbasin))
       do n = 1, model%ocean_data%nbasin
          model%ocean_data%basin(n) = n
       enddo
       allocate(model%scalars%iarea_basin(model%ocean_data%nbasin))
       allocate(model%scalars%iareag_basin(model%ocean_data%nbasin))
       allocate(model%scalars%iareaf_basin(model%ocean_data%nbasin))
       allocate(model%scalars%ivol_basin(model%ocean_data%nbasin))
       allocate(model%scalars%ivol_above_flotation_basin(model%ocean_data%nbasin))
       allocate(model%scalars%imass_basin(model%ocean_data%nbasin))
       allocate(model%scalars%imass_above_flotation_basin(model%ocean_data%nbasin))
       allocate(model%scalars%icap_area_basin(model%ocean_data%nbasin))
       allocate(model%scalars%icap_vol_basin(model%ocean_data%nbasin))
       allocate(model%scalars%bmlt_float_basin(model%ocean_data%nbasin))
       allocate(model%scalars%bmlt_float_target_basin(model%ocean_data%nbasin))
    endif

    ! climate arrays
    call coordsystem_allocate(model%grid%ice_grid, model%climate%acab)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%acab_applied)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%acab_applied_tavg)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%acab_cpl)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%acab_cpl_tavg)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%acab_unapplied_glad)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%acab_unapplied_glad_tavg)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%acab_unapplied)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%acab_unapplied_tavg)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%smb)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%smb_anomaly)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%smb_corrected)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%smb_obs)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%artm)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%artm_anomaly)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%artm_corrected)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%snow)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%snow_anomaly)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%snow_corrected)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%precip)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%precip_anomaly)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%precip_corrected)
    call coordsystem_allocate(model%grid%ice_grid, model%climate%overwrite_acab_mask)

    if (model%climate%smb_input_function == SMB_INPUT_FUNCTION_XY_GRADZ) then
       call coordsystem_allocate(model%grid%ice_grid, model%climate%smb_ref)
       call coordsystem_allocate(model%grid%ice_grid, model%climate%smb_gradz)
       if (.not.associated(model%climate%usrf_ref)) &
            call coordsystem_allocate(model%grid%ice_grid, model%climate%usrf_ref)
    elseif (model%climate%smb_input_function == SMB_INPUT_FUNCTION_XYZ) then
       call coordsystem_allocate(model%grid%ice_grid, model%climate%nzatm, model%climate%smb_3d)
       if (.not.associated(model%climate%zatm)) allocate(model%climate%zatm(model%climate%nzatm))
    elseif (model%climate%smb_input_function == SMB_INPUT_FUNCTION_PDD) then
       call coordsystem_allocate(model%grid%ice_grid, model%climate%ablation)
    endif

    ! Note: Typically, smb_input_function and artm_input_function will have the same value.
    !       If both use a lapse rate, they will share the array usrf_ref
    !       If both are 3d, they will share the array smb_levels.
    if (model%climate%artm_input_function == ARTM_INPUT_FUNCTION_XY_GRADZ) then
       call coordsystem_allocate(model%grid%ice_grid, model%climate%artm_ref)
       call coordsystem_allocate(model%grid%ice_grid, model%climate%artm_gradz)
       if (.not.associated(model%climate%usrf_ref)) then
          call coordsystem_allocate(model%grid%ice_grid, model%climate%usrf_ref)
       endif
    elseif (model%climate%smb_input_function == ARTM_INPUT_FUNCTION_XYZ) then
       call coordsystem_allocate(model%grid%ice_grid, model%climate%nzatm, model%climate%artm_3d)
       if (.not.associated(model%climate%zatm)) allocate(model%climate%zatm(model%climate%nzatm))
    elseif (model%climate%artm_input_function == ARTM_INPUT_FUNCTION_XY_LAPSE) then
       call coordsystem_allocate(model%grid%ice_grid, model%climate%artm_ref)
       if (.not.associated(model%climate%usrf_ref)) then
          call coordsystem_allocate(model%grid%ice_grid, model%climate%usrf_ref)
       endif
    endif

    ! calving arrays
    call coordsystem_allocate(model%grid%ice_grid, model%calving%calving_thck)
    call coordsystem_allocate(model%grid%ice_grid, model%calving%calving_rate)
    call coordsystem_allocate(model%grid%ice_grid, model%calving%calving_rate_tavg)
    if (.not. model%calving%subgrid_calving_front) then
       call coordsystem_allocate(model%grid%ice_grid, model%calving%calving_mask)
    else
       call coordsystem_allocate(model%grid%ice_grid, model%calving%subgrid_calving_mask)
    endif
    call coordsystem_allocate(model%grid%ice_grid, model%calving%calving_front_mask)
    call coordsystem_allocate(model%grid%ice_grid, model%calving%beyond_cf_mask)
    call coordsystem_allocate(model%grid%ice_grid, model%calving%thck_effective)
    call coordsystem_allocate(model%grid%ice_grid, model%calving%effective_areafrac)
    call coordsystem_allocate(model%grid%ice_grid, model%calving%tau_eigen1)
    call coordsystem_allocate(model%grid%ice_grid, model%calving%tau_eigen2)
    call coordsystem_allocate(model%grid%ice_grid, model%calving%eps_eigen1)
    call coordsystem_allocate(model%grid%ice_grid, model%calving%eps_eigen2)
    if (model%calving%calving_law == CALVING_LAW_DAMAGE) then
       call coordsystem_allocate(model%grid%ice_grid, upn-1, model%calving%damage)
    else
       ! allocate with size 1, since they need to be allocated to be passed to calving subroutine
       allocate(model%calving%damage(1,1,1))
    endif
    if (model%calving%calvingmip_domain /= CALVINGMIP_DOMAIN_NONE) then
       model%calving%naxis = 8
       allocate(model%calving%axis(model%calving%naxis))
       do n = 1, model%calving%naxis
          model%calving%axis(n) = n
       enddo
       allocate(model%calving%cf_locx(model%calving%naxis))
       allocate(model%calving%cf_locy(model%calving%naxis))
       allocate(model%calving%cf_radius(model%calving%naxis))
       allocate(model%calving%cf_thck(model%calving%naxis))
       allocate(model%calving%cf_uvel(model%calving%naxis))
       allocate(model%calving%cf_vvel(model%calving%naxis))
    endif

    ! lateral melt arrays
    ! Note: latmelt_thck and latmelt_rate are used in mass balance diagnostics
    !       even when lateral melt is turned off
    call coordsystem_allocate(model%grid%ice_grid, model%lateral_melt%latmelt_thck)
    call coordsystem_allocate(model%grid%ice_grid, model%lateral_melt%latmelt_rate)
    call coordsystem_allocate(model%grid%ice_grid, model%lateral_melt%latmelt_rate_tavg)
    if (model%lateral_melt%lateral_melt_method /= LATERAL_MELT_METHOD_NONE) then
       call coordsystem_allocate(model%grid%ice_grid, model%lateral_melt%runoff_applied)
       call coordsystem_allocate(model%grid%ice_grid, model%lateral_melt%thermal_forcing_applied)
       call coordsystem_allocate(model%grid%ice_grid, model%lateral_melt%melt_front_mask)
    endif

    ! isostasy arrays

    ! Note: relx needs to be allocated always
    call coordsystem_allocate(model%grid%ice_grid, model%isostasy%relx)
    if (model%isostasy%compute_isostasy) then
       call coordsystem_allocate(model%grid%ice_grid, model%isostasy%load)
       call coordsystem_allocate(model%grid%ice_grid, model%isostasy%load_factors)
       call coordsystem_allocate(model%grid%ice_grid, model%isostasy%tau_relax)
    endif

  end subroutine glissade_allocarr

!=====================================================================
  
  subroutine glissade_deallocarr(model)

    !> deallocate model arrays
    !TODO - Check that all arrays allocated above are deallocated here.

    implicit none
    type(glissade_global_type),intent(inout) :: model

    ! horizontal coordinates

    if (associated(model%grid%x0)) &
        deallocate(model%grid%x0) 
    if (associated(model%grid%y0)) &
        deallocate(model%grid%y0) 
    if (associated(model%grid%x1)) &
        deallocate(model%grid%x1) 
    if (associated(model%grid%y1)) &
        deallocate(model%grid%y1) 
    if (associated(model%grid%x0_global)) &
        deallocate(model%grid%x0_global)
    if (associated(model%grid%y0_global)) &
        deallocate(model%grid%y0_global)
    if (associated(model%grid%x1_global)) &
        deallocate(model%grid%x1_global)
    if (associated(model%grid%y1_global)) &
        deallocate(model%grid%y1_global)
    if (associated(model%grid%lat)) &
        deallocate(model%grid%lat)
    if (associated(model%grid%lon)) &
        deallocate(model%grid%lon)
    if (associated(model%grid%ice_domain_mask)) &
        deallocate(model%grid%ice_domain_mask)
    if (associated(model%grid%global_edge_mask)) &
        deallocate(model%grid%global_edge_mask)

    ! vertical sigma coordinates

    if (associated(model%grid%sigma)) &
        deallocate(model%grid%sigma)
    if (associated(model%grid%stagsigma)) &
        deallocate(model%grid%stagsigma)
    if (associated(model%grid%stagwbndsigma)) &
        deallocate(model%grid%stagwbndsigma)

    ! temperature arrays

    if (associated(model%temperature%temp)) &
        deallocate(model%temperature%temp)
    if (associated(model%temperature%tempunstag)) &
        deallocate(model%temperature%tempunstag)
    if (associated(model%temperature%bheatflx)) &
        deallocate(model%temperature%bheatflx)
    if (associated(model%temperature%bpmp)) &
        deallocate(model%temperature%bpmp)
    if (associated(model%temperature%btemp)) &
        deallocate(model%temperature%btemp)
    if (associated(model%temperature%btemp_ground)) &
        deallocate(model%temperature%btemp_ground)
    if (associated(model%temperature%btemp_float)) &
        deallocate(model%temperature%btemp_float)
    if (associated(model%temperature%bfricflx)) &
        deallocate(model%temperature%bfricflx)
    if (associated(model%temperature%bhydroflx)) &
        deallocate(model%temperature%bhydroflx)
    if (associated(model%temperature%ucondflx)) &
        deallocate(model%temperature%ucondflx)
    if (associated(model%temperature%waterfrac)) &
        deallocate(model%temperature%waterfrac)
    if (associated(model%temperature%enthalpy)) &
        deallocate(model%temperature%enthalpy)
    if (associated(model%temperature%flwa)) &
        deallocate(model%temperature%flwa)
    if (associated(model%temperature%dissip)) &
        deallocate(model%temperature%dissip)
    if (associated(model%temperature%flow_enhancement_factor)) &
        deallocate(model%temperature%flow_enhancement_factor)

    ! velocity arrays

    if (associated(model%velocity%uvel)) &
        deallocate(model%velocity%uvel)
    if (associated(model%velocity%vvel)) &
        deallocate(model%velocity%vvel)
    if (associated(model%velocity%uvel_2d)) &
        deallocate(model%velocity%uvel_2d)
    if (associated(model%velocity%vvel_2d)) &
        deallocate(model%velocity%vvel_2d)
    if (associated(model%velocity%velnorm)) &
        deallocate(model%velocity%velnorm)
    if (associated(model%velocity%wvel)) &
        deallocate(model%velocity%wvel)
    if (associated(model%velocity%btrc)) &
        deallocate(model%velocity%btrc)
    if (associated(model%velocity%resid_u)) &
        deallocate(model%velocity%resid_u)
    if (associated(model%velocity%resid_v)) &
        deallocate(model%velocity%resid_v)
    if (associated(model%velocity%rhs_u)) &
        deallocate(model%velocity%rhs_u)
    if (associated(model%velocity%rhs_v)) &
        deallocate(model%velocity%rhs_v)
    if (associated(model%velocity%uvel_extend)) &
        deallocate(model%velocity%uvel_extend)
    if (associated(model%velocity%vvel_extend)) &
        deallocate(model%velocity%vvel_extend)
    if (associated(model%velocity%uvel_2d_extend)) &
        deallocate(model%velocity%uvel_2d_extend)
    if (associated(model%velocity%vvel_2d_extend)) &
        deallocate(model%velocity%vvel_2d_extend)
    if (associated(model%velocity%uvel_mean)) &
        deallocate(model%velocity%uvel_mean)
    if (associated(model%velocity%vvel_mean)) &
        deallocate(model%velocity%vvel_mean)
    if (associated(model%velocity%velnorm_mean)) &
        deallocate(model%velocity%velnorm_mean)
    if (associated(model%velocity%usfc_obs)) &
        deallocate(model%velocity%usfc_obs)
    if (associated(model%velocity%vsfc_obs)) &
        deallocate(model%velocity%vsfc_obs)
    if (associated(model%velocity%velo_sfc_obs)) &
        deallocate(model%velocity%velo_sfc_obs)
    if (associated(model%velocity%velo_sfc)) &
        deallocate(model%velocity%velo_sfc)

!!    if (associated(model%velocity%velmask)) &   ! no longer used
!!       deallocate(model%velocity%velmask) 

    if (associated(model%velocity%beta)) &
        deallocate(model%velocity%beta)
    if (associated(model%velocity%unstagbeta)) &
        deallocate(model%velocity%unstagbeta)
    if (associated(model%velocity%beta_internal)) &
        deallocate(model%velocity%beta_internal)
    if (associated(model%velocity%kinbcmask)) &
        deallocate(model%velocity%kinbcmask)
    if (associated(model%velocity%umask_no_penetration)) &
        deallocate(model%velocity%umask_no_penetration)
    if (associated(model%velocity%vmask_no_penetration)) &
        deallocate(model%velocity%vmask_no_penetration)
    if (associated(model%velocity%divu)) &
        deallocate(model%velocity%divu)
    if (associated(model%velocity%shear)) &
        deallocate(model%velocity%shear)

    ! higher-order stress arrays

    if (associated(model%stress%efvs)) &
        deallocate(model%stress%efvs)
    if (associated(model%stress%efvs_vertavg)) &
        deallocate(model%stress%efvs_vertavg)
    if (associated(model%stress%tau%scalar)) &
        deallocate(model%stress%tau%scalar)
    if (associated(model%stress%tau%xz)) &
        deallocate(model%stress%tau%xz)
    if (associated(model%stress%tau%yz)) &
        deallocate(model%stress%tau%yz)
    if (associated(model%stress%tau%xx)) &
        deallocate(model%stress%tau%xx)
    if (associated(model%stress%tau%yy)) &
        deallocate(model%stress%tau%yy)
    if (associated(model%stress%tau%xy)) &
        deallocate(model%stress%tau%xy)
    if (associated(model%velocity%strain_rate%scalar)) &
        deallocate(model%velocity%strain_rate%scalar)
    if (associated(model%velocity%strain_rate%xz)) &
        deallocate(model%velocity%strain_rate%xz)
    if (associated(model%velocity%strain_rate%yz)) &
        deallocate(model%velocity%strain_rate%yz)
    if (associated(model%velocity%strain_rate%xx)) &
        deallocate(model%velocity%strain_rate%xx)
    if (associated(model%velocity%strain_rate%yy)) &
        deallocate(model%velocity%strain_rate%yy)
    if (associated(model%velocity%strain_rate%xy)) &
        deallocate(model%velocity%strain_rate%xy)
    if (associated(model%stress%btractx)) &
        deallocate(model%stress%btractx)
    if (associated(model%stress%btracty)) &
        deallocate(model%stress%btracty)
    if (associated(model%stress%btract)) &
        deallocate(model%stress%btract)
    if (associated(model%stress%btractx_extend)) &
        deallocate(model%stress%btractx_extend)
    if (associated(model%stress%btracty_extend)) &
        deallocate(model%stress%btracty_extend)
    if (associated(model%stress%taudx)) &
        deallocate(model%stress%taudx)
    if (associated(model%stress%taudy)) &
        deallocate(model%stress%taudy)

    ! basal hydrology arrays
    if (associated(model%basal_hydro%bwat)) &
        deallocate(model%basal_hydro%bwat)
    if (associated(model%basal_hydro%bwatflx)) &
        deallocate(model%basal_hydro%bwatflx)
    if (associated(model%basal_hydro%bmlt_hydro)) &
        deallocate(model%basal_hydro%bmlt_hydro)
    if (associated(model%basal_hydro%bwat_diag)) &
        deallocate(model%basal_hydro%bwat_diag)
    if (associated(model%basal_hydro%head)) &
        deallocate(model%basal_hydro%head)
    if (associated(model%basal_hydro%grad_head)) &
        deallocate(model%basal_hydro%grad_head)

    ! basal physics arrays
    if (associated(model%basal_physics%bpmp_mask)) &
        deallocate(model%basal_physics%bpmp_mask)
    if (associated(model%basal_physics%effecpress)) &
        deallocate(model%basal_physics%effecpress)
    if (associated(model%basal_physics%effecpress_stag)) &
        deallocate(model%basal_physics%effecpress_stag)
    if (associated(model%basal_physics%c_space_factor)) &
        deallocate(model%basal_physics%c_space_factor)
    if (associated(model%basal_physics%c_space_factor_stag)) &
        deallocate(model%basal_physics%c_space_factor_stag)
    if (associated(model%basal_physics%mintauf)) &
       deallocate(model%basal_physics%mintauf)

    ! basal melt arrays

    if (associated(model%basal_melt%bmlt)) &
        deallocate(model%basal_melt%bmlt)
    if (associated(model%basal_melt%bmlt_applied)) &
        deallocate(model%basal_melt%bmlt_applied)
    if (associated(model%basal_melt%bmlt_applied_tavg)) &
        deallocate(model%basal_melt%bmlt_applied_tavg)
    if (associated(model%basal_melt%bmlt_ground)) &
        deallocate(model%basal_melt%bmlt_ground)
    if (associated(model%basal_melt%bmlt_float)) &
        deallocate(model%basal_melt%bmlt_float)
    if (associated(model%basal_melt%bmlt_float_target)) &
        deallocate(model%basal_melt%bmlt_float_target)
    if (associated(model%basal_melt%bmlt_float_external)) &
        deallocate(model%basal_melt%bmlt_float_external)
    if (associated(model%basal_melt%bmlt_float_anomaly)) &
        deallocate(model%basal_melt%bmlt_float_anomaly)
    if (associated(model%basal_melt%warm_ocean_mask)) &
        deallocate(model%basal_melt%warm_ocean_mask)
    if (associated(model%basal_melt%thermal_forcing_mask)) &
        deallocate(model%basal_melt%thermal_forcing_mask)

    ! ocean data arrays
    if (associated(model%ocean_data%basin_number)) &
        deallocate(model%ocean_data%basin_number)
    if (associated(model%ocean_data%deltaT_ocn)) &
        deallocate(model%ocean_data%deltaT_ocn)
    if (associated(model%ocean_data%thermal_forcing)) &
        deallocate(model%ocean_data%thermal_forcing)
    if (associated(model%ocean_data%thermal_forcing_lsrf)) &
        deallocate(model%ocean_data%thermal_forcing_lsrf)
    if (associated(model%ocean_data%bmb_float)) &
        deallocate(model%ocean_data%bmb_float)

    ! glacier arrays
    if (associated(model%glaciers%glacierid)) &
        deallocate(model%glaciers%glacierid)
    if (associated(model%glaciers%rgi_glacier_id)) &
        deallocate(model%glaciers%rgi_glacier_id)
    if (associated(model%glaciers%cism_glacier_id)) &
        deallocate(model%glaciers%cism_glacier_id)
    if (associated(model%glaciers%cism_glacier_id_init)) &
        deallocate(model%glaciers%cism_glacier_id_init)
    if (associated(model%glaciers%smb_glacier_id)) &
        deallocate(model%glaciers%smb_glacier_id)
    if (associated(model%glaciers%smb_glacier_id_init)) &
        deallocate(model%glaciers%smb_glacier_id_init)
    if (associated(model%glaciers%cism_to_rgi_glacier_id)) &
        deallocate(model%glaciers%cism_to_rgi_glacier_id)
    if (associated(model%glaciers%dthck_dt_annmean)) &
        deallocate(model%glaciers%dthck_dt_annmean)
    if (associated(model%glaciers%snow_annmean)) &
        deallocate(model%glaciers%snow_annmean)
    if (associated(model%glaciers%Tpos_annmean)) &
        deallocate(model%glaciers%Tpos_annmean)
    if (associated(model%glaciers%snow_rgi_annmean)) &
        deallocate(model%glaciers%snow_rgi_annmean)
    if (associated(model%glaciers%Tpos_rgi_annmean)) &
        deallocate(model%glaciers%Tpos_rgi_annmean)
    if (associated(model%glaciers%snow_recent_annmean)) &
        deallocate(model%glaciers%snow_recent_annmean)
    if (associated(model%glaciers%Tpos_recent_annmean)) &
        deallocate(model%glaciers%Tpos_recent_annmean)
    if (associated(model%glaciers%smb_applied_annmean)) &
        deallocate(model%glaciers%smb_applied_annmean)
    if (associated(model%glaciers%smb_obs)) &
        deallocate(model%glaciers%smb_obs)
    if (associated(model%glaciers%area)) &
        deallocate(model%glaciers%area)
    if (associated(model%glaciers%volume)) &
        deallocate(model%glaciers%volume)
    if (associated(model%glaciers%area_init)) &
        deallocate(model%glaciers%area_init)
    if (associated(model%glaciers%volume_init)) &
        deallocate(model%glaciers%volume_init)
    if (associated(model%glaciers%area_init_extent)) &
        deallocate(model%glaciers%area_init_extent)
    if (associated(model%glaciers%volume_init_extent)) &
        deallocate(model%glaciers%volume_init_extent)
    if (associated(model%glaciers%area_target)) &
        deallocate(model%glaciers%area_target)
    if (associated(model%glaciers%volume_target)) &
        deallocate(model%glaciers%volume_target)
    if (associated(model%glaciers%mu_star)) &
        deallocate(model%glaciers%mu_star)
    if (associated(model%glaciers%alpha_snow)) &
        deallocate(model%glaciers%alpha_snow)
    if (associated(model%glaciers%beta_artm)) &
        deallocate(model%glaciers%beta_artm)
    if (associated(model%glaciers%smb)) &
        deallocate(model%glaciers%smb)
    if (associated(model%glaciers%usrf_target)) &
        deallocate(model%glaciers%usrf_target)
    if (associated(model%glaciers%thck_target)) &
        deallocate(model%glaciers%thck_target)
    if (associated(model%glaciers%smb_rgi)) &
        deallocate(model%glaciers%smb_rgi)
    if (associated(model%glaciers%delta_usrf_rgi)) &
        deallocate(model%glaciers%delta_usrf_rgi)
    if (associated(model%glaciers%smb_recent)) &
        deallocate(model%glaciers%smb_recent)
    if (associated(model%glaciers%delta_usrf_recent)) &
        deallocate(model%glaciers%delta_usrf_recent)
    if (associated(model%glaciers%boundary_mask)) &
        deallocate(model%glaciers%boundary_mask)

    ! inversion arrays
    if (associated(model%basal_physics%powerlaw_c)) &
        deallocate(model%basal_physics%powerlaw_c)
    if (associated(model%basal_physics%coulomb_c)) &
        deallocate(model%basal_physics%coulomb_c)
    if (associated(model%inversion%floating_thck_target)) &
        deallocate(model%inversion%floating_thck_target)
    if (associated(model%inversion%grounded_thck_target)) &
        deallocate(model%inversion%grounded_thck_target)

    ! geometry arrays

    if (associated(model%geometry%thck)) &
        deallocate(model%geometry%thck)
    if (associated(model%geometry%stagthck)) &
        deallocate(model%geometry%stagthck)
    if (associated(model%geometry%usrf)) &
        deallocate(model%geometry%usrf)
    if (associated(model%geometry%lsrf)) &
        deallocate(model%geometry%lsrf)
    if (associated(model%geometry%topg)) &
        deallocate(model%geometry%topg)
    if (associated(model%geometry%usrf_obs)) &
        deallocate(model%geometry%usrf_obs)
    if (associated(model%geometry%dthck_dt)) &
        deallocate(model%geometry%dthck_dt)
    if (associated(model%geometry%dthck_dt_obs)) &
        deallocate(model%geometry%dthck_dt_obs)
    if (associated(model%geometry%dthck_dt_obs_basin)) &
        deallocate(model%geometry%dthck_dt_obs_basin)
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

    if (associated(model%geometry%tracers)) &
        deallocate(model%geometry%tracers)
    if (associated(model%geometry%tracers_usrf)) &
        deallocate(model%geometry%tracers_usrf)
    if (associated(model%geometry%tracers_lsrf)) &
        deallocate(model%geometry%tracers_lsrf)

    if (associated(model%mass_flux%sfc_mbal_flux)) &
        deallocate(model%mass_flux%sfc_mbal_flux)
    if (associated(model%mass_flux%sfc_mbal_flux_tavg)) &
        deallocate(model%mass_flux%sfc_mbal_flux_tavg)
    if (associated(model%mass_flux%basal_mbal_flux)) &
        deallocate(model%mass_flux%basal_mbal_flux)
    if (associated(model%mass_flux%basal_mbal_flux_tavg)) &
        deallocate(model%mass_flux%basal_mbal_flux_tavg)
    if (associated(model%mass_flux%calving_flux)) &
        deallocate(model%mass_flux%calving_flux)
    if (associated(model%mass_flux%calving_flux_tavg)) &
        deallocate(model%mass_flux%calving_flux_tavg)
    if (associated(model%mass_flux%latmelt_flux)) &
        deallocate(model%mass_flux%latmelt_flux)
    if (associated(model%mass_flux%latmelt_flux_tavg)) &
        deallocate(model%mass_flux%latmelt_flux_tavg)
    if (associated(model%mass_flux%removal_flux)) &
        deallocate(model%mass_flux%removal_flux)
    if (associated(model%mass_flux%removal_flux_tavg)) &
        deallocate(model%mass_flux%removal_flux_tavg)
    if (associated(model%mass_flux%gl_flux_east)) &
        deallocate(model%mass_flux%gl_flux_east)
    if (associated(model%mass_flux%gl_flux_north)) &
        deallocate(model%mass_flux%gl_flux_north)
    if (associated(model%mass_flux%gl_flux)) &
        deallocate(model%mass_flux%gl_flux)
    if (associated(model%mass_flux%gl_flux_tavg)) &
        deallocate(model%mass_flux%gl_flux_tavg)

    if (associated(model%geometry%ice_mask)) &
       deallocate(model%geometry%ice_mask)
    if (associated(model%geometry%ice_mask_stag)) &
       deallocate(model%geometry%ice_mask_stag)
    if (associated(model%geometry%floating_mask)) &
       deallocate(model%geometry%floating_mask)
    if (associated(model%geometry%grounded_mask)) &
       deallocate(model%geometry%grounded_mask)
    if (associated(model%geometry%marine_connection_mask)) &
       deallocate(model%geometry%marine_connection_mask)
    if (associated(model%geometry%ice_sheet_mask)) &
       deallocate(model%geometry%ice_sheet_mask)
    if (associated(model%geometry%ice_cap_mask)) &
       deallocate(model%geometry%ice_cap_mask)
    if (associated(model%geometry%ice_fraction_retreat_mask)) &
       deallocate(model%geometry%ice_fraction_retreat_mask)
    if (associated(model%geometry%removal_thck)) &
       deallocate(model%geometry%removal_thck)
    if (associated(model%geometry%removal_rate)) &
       deallocate(model%geometry%removal_rate)
    if (associated(model%geometry%removal_rate_tavg)) &
       deallocate(model%geometry%removal_rate_tavg)
    if (associated(model%geometry%reference_thck)) &
       deallocate(model%geometry%reference_thck)

    if (associated(model%scalars%iarea_basin)) &
         deallocate(model%scalars%iarea_basin)
    if (associated(model%scalars%iareag_basin)) &
         deallocate(model%scalars%iareag_basin)
    if (associated(model%scalars%iareaf_basin)) &
         deallocate(model%scalars%iareaf_basin)
    if (associated(model%scalars%ivol_basin)) &
         deallocate(model%scalars%ivol_basin)
    if (associated(model%scalars%ivol_above_flotation_basin)) &
         deallocate(model%scalars%ivol_above_flotation_basin)
    if (associated(model%scalars%imass_basin)) &
         deallocate(model%scalars%imass_basin)
    if (associated(model%scalars%imass_above_flotation_basin)) &
         deallocate(model%scalars%imass_above_flotation_basin)
    if (associated(model%scalars%icap_area_basin)) &
         deallocate(model%scalars%icap_area_basin)
    if (associated(model%scalars%icap_vol_basin)) &
         deallocate(model%scalars%icap_vol_basin)
    if (associated(model%scalars%bmlt_float_basin)) &
         deallocate(model%scalars%bmlt_float_basin)
    if (associated(model%scalars%bmlt_float_target_basin)) &
         deallocate(model%scalars%bmlt_float_target_basin)

    if (associated(model%geometry%ice_age)) &
        deallocate(model%geometry%ice_age)
    if (associated(model%geometry%thck_old)) &
        deallocate(model%geometry%thck_old)
    if (associated(model%geometry%tracers)) &
        deallocate(model%geometry%tracers)
    if (associated(model%geometry%f_flotation)) &
        deallocate(model%geometry%f_flotation)
    if (associated(model%geometry%f_ground)) &
        deallocate(model%geometry%f_ground)
    if (associated(model%geometry%f_ground_cell)) &
        deallocate(model%geometry%f_ground_cell)

    ! climate arrays

    if (associated(model%climate%acab)) &
        deallocate(model%climate%acab)
    if (associated(model%climate%acab_applied)) &
        deallocate(model%climate%acab_applied)
    if (associated(model%climate%acab_applied_tavg)) &
        deallocate(model%climate%acab_applied_tavg)
    if (associated(model%climate%acab_cpl)) &
        deallocate(model%climate%acab_cpl)
    if (associated(model%climate%acab_cpl_tavg)) &
        deallocate(model%climate%acab_cpl_tavg)
    if (associated(model%climate%acab_unapplied_glad)) &
        deallocate(model%climate%acab_unapplied_glad)
    if (associated(model%climate%acab_unapplied_glad_tavg)) &
        deallocate(model%climate%acab_unapplied_glad_tavg)
    if (associated(model%climate%acab_unapplied)) &
        deallocate(model%climate%acab_unapplied)
    if (associated(model%climate%acab_unapplied_tavg)) &
        deallocate(model%climate%acab_unapplied_tavg)
    if (associated(model%climate%smb)) &
        deallocate(model%climate%smb)
    if (associated(model%climate%smb_anomaly)) &
        deallocate(model%climate%smb_anomaly)
    if (associated(model%climate%smb_corrected)) &
        deallocate(model%climate%smb_corrected)
    if (associated(model%climate%artm)) &
        deallocate(model%climate%artm)
    if (associated(model%climate%artm_anomaly)) &
        deallocate(model%climate%artm_anomaly)
    if (associated(model%climate%artm_corrected)) &
        deallocate(model%climate%artm_corrected)
    if (associated(model%climate%snow)) &
        deallocate(model%climate%snow)
    if (associated(model%climate%snow_anomaly)) &
        deallocate(model%climate%snow_anomaly)
    if (associated(model%climate%snow_corrected)) &
        deallocate(model%climate%snow_corrected)
    if (associated(model%climate%precip)) &
        deallocate(model%climate%precip)
    if (associated(model%climate%precip_anomaly)) &
        deallocate(model%climate%precip_anomaly)
    if (associated(model%climate%precip_corrected)) &
        deallocate(model%climate%precip_corrected)
    if (associated(model%climate%overwrite_acab_mask)) &
        deallocate(model%climate%overwrite_acab_mask)
    if (associated(model%climate%smb_ref)) &
        deallocate(model%climate%smb_ref)
    if (associated(model%climate%smb_gradz)) &
        deallocate(model%climate%smb_gradz)
    if (associated(model%climate%usrf_ref)) &
        deallocate(model%climate%usrf_ref)
    if (associated(model%climate%artm_ref)) &
        deallocate(model%climate%artm_ref)
    if (associated(model%climate%artm_gradz)) &
        deallocate(model%climate%artm_gradz)
    if (associated(model%climate%smb_3d)) &
        deallocate(model%climate%smb_3d)
    if (associated(model%climate%artm_3d)) &
        deallocate(model%climate%artm_3d)
    if (associated(model%climate%smb_obs)) &
        deallocate(model%climate%smb_obs)
    if (associated(model%climate%ablation)) &
        deallocate(model%climate%ablation)

    ! calving arrays
    if (associated(model%calving%calving_thck)) &
        deallocate(model%calving%calving_thck)
    if (associated(model%calving%calving_rate)) &
        deallocate(model%calving%calving_rate)
    if (associated(model%calving%calving_rate_tavg)) &
        deallocate(model%calving%calving_rate_tavg)
    if (associated(model%calving%calving_mask)) &
        deallocate(model%calving%calving_mask)
    if (associated(model%calving%subgrid_calving_mask)) &
        deallocate(model%calving%subgrid_calving_mask)
    if (associated(model%calving%calving_front_mask)) &
        deallocate(model%calving%calving_front_mask)
    if (associated(model%calving%beyond_cf_mask)) &
        deallocate(model%calving%beyond_cf_mask)
    if (associated(model%calving%thck_effective)) &
        deallocate(model%calving%thck_effective)
    if (associated(model%calving%effective_areafrac)) &
        deallocate(model%calving%effective_areafrac)
    if (associated(model%calving%tau_eigen1)) &
        deallocate(model%calving%tau_eigen1)
    if (associated(model%calving%tau_eigen2)) &
        deallocate(model%calving%tau_eigen2)
    if (associated(model%calving%eps_eigen1)) &
        deallocate(model%calving%eps_eigen1)
    if (associated(model%calving%eps_eigen2)) &
        deallocate(model%calving%eps_eigen2)
    if (associated(model%calving%damage)) &
        deallocate(model%calving%damage)
    if (associated(model%calving%axis)) &
        deallocate(model%calving%axis)
    if (associated(model%calving%cf_locx)) &
        deallocate(model%calving%cf_locx)
    if (associated(model%calving%cf_locy)) &
        deallocate(model%calving%cf_locy)
    if (associated(model%calving%cf_radius)) &
        deallocate(model%calving%cf_radius)
    if (associated(model%calving%cf_thck)) &
        deallocate(model%calving%cf_thck)
    if (associated(model%calving%cf_uvel)) &
        deallocate(model%calving%cf_uvel)
    if (associated(model%calving%cf_vvel)) &
        deallocate(model%calving%cf_vvel)

    ! lateral melt arrays
    if (associated(model%lateral_melt%runoff_applied)) &
        deallocate(model%lateral_melt%runoff_applied)
    if (associated(model%lateral_melt%thermal_forcing_applied)) &
        deallocate(model%lateral_melt%thermal_forcing_applied)
    if (associated(model%lateral_melt%latmelt_thck)) &
        deallocate(model%lateral_melt%latmelt_thck)
    if (associated(model%lateral_melt%latmelt_rate)) &
        deallocate(model%lateral_melt%latmelt_rate)
    if (associated(model%lateral_melt%latmelt_rate_tavg)) &
        deallocate(model%lateral_melt%latmelt_rate_tavg)
    if (associated(model%lateral_melt%melt_front_mask)) &
        deallocate(model%lateral_melt%melt_front_mask)

    ! isostasy arrays

    if (associated(model%isostasy%relx)) &
        deallocate(model%isostasy%relx)
    if (associated(model%isostasy%load)) &
        deallocate(model%isostasy%load)
    if (associated(model%isostasy%load_factors)) &
        deallocate(model%isostasy%load_factors)
    if (associated(model%isostasy%tau_relax)) &
        deallocate(model%isostasy%tau_relax)

    ! projection arrays
    if (associated(model%projection%stere)) then
       if (associated(model%projection%stere%area_factor)) then
          deallocate(model%projection%stere%area_factor)
       endif
    endif

  end subroutine glissade_deallocarr

!=====================================================================
  
  ! some accessor functions
  function get_dew(model)
    !> return scaled x node spacing
    implicit none
    real(dp) :: get_dew
    type(glissade_global_type) :: model

    get_dew = model%grid%dew
  end function get_dew

  function get_dns(model)
    !> return scaled y node spacing
    implicit none
    real(dp) :: get_dns
    type(glissade_global_type) :: model

    get_dns = model%grid%dns
  end function get_dns

  function get_tstart(model)
    !> return start time
    implicit none
    real(dp) :: get_tstart
    type(glissade_global_type) :: model
    
    get_tstart = model%time%tstart
  end function get_tstart

  function get_tend(model)
    !> return end time
    implicit none
    real(dp) :: get_tend
    type(glissade_global_type) :: model
    
    get_tend = model%time%tend
  end function get_tend

  function get_tinc(model)
    !> return time increment
    implicit none
    real(dp) :: get_tinc
    type(glissade_global_type) :: model
    
    get_tinc = model%time%tinc
  end function get_tinc

  function get_ewn(model)
    !> get number of nodes in x dir
    implicit none
    integer get_ewn
    type(glissade_global_type) :: model

    get_ewn = model%grid%ewn
  end function get_ewn

  function get_nsn(model)
    !> get number of nodes in y dir
    implicit none
    integer get_nsn
    type(glissade_global_type) :: model

    get_nsn = model%grid%nsn
  end function get_nsn

  function get_nzocn(model)
    !> get number of ocean layers
    implicit none
    integer get_nzocn
    type(glissade_global_type) :: model

    get_nzocn = model%ocean_data%nzocn
  end function get_nzocn

  function get_nzatm(model)
    !> get number of atmosphere layers
    implicit none
    integer get_nzatm
    type(glissade_global_type) :: model

    get_nzatm = model%climate%nzatm
  end function get_nzatm

  subroutine set_time(model,time)
    !> Set the model time counter --- useful for
    !> fractional year output
    implicit none
    type(glissade_global_type) :: model
    real(dp) :: time

    model%time%time = time

  end subroutine set_time

!=====================================================================  

end module glissade_types

!=====================================================================
