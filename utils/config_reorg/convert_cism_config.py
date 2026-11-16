#!/usr/bin/env python3
"""Convert CISM (Glissade) config files, or CESM user_nl_cism files, to the reorganized config format.

The config reorganization moves keys out of the long [options], [ho_options] and [parameters]
sections into short sections by topic (e.g. [calving], [climate], [basal_physics]), and renames
some keys (e.g. marine_margin -> calving_law).  Option values are not changed.

Usage:
    convert_cism_config.py FILE [FILE ...]              write FILE.converted next to each FILE
    convert_cism_config.py FILE -o NEWFILE              write one converted file to NEWFILE
    convert_cism_config.py FILE [FILE ...] --in-place   overwrite each FILE
    convert_cism_config.py FILE [FILE ...] --check      report what would change; exit 1 if anything would
    convert_cism_config.py --self-test                  run the built-in tests

A file with no [section] lines is treated as a CESM user_nl_cism file, in which only key names change
(the CESM namelist has no sections).  Use --user-nl to force this.

How a config file is converted:
  - Each key listed in MOVES is moved to its new section, under its new name, with its value and
    inline comment.  Comment lines directly above the key (with no blank line in between) move with it.
  - A new section is created if needed.
  - An old section left with no keys is removed, together with any comment lines left in it.
  - The Glissade sections are put in the standard order (SECTION_ORDER: the order of the log file),
    followed by the [CF ...] sections and any other sections, in their original order.  Comment lines
    above a section header (after a blank line) move with that section.  Section contents are not
    changed.
  - Keys listed in REMOVED (keys the model no longer reads) are deleted; the report names the replacement,
    if there is one, e.g. 'deleted [isostasy] lithosphere_period: removed; see load_update_interval'.
  - Other keys are left where they are; the model warns about keys it does not read.
  - Converting a file that is already converted changes nothing.
  - Glide configs ([options] dycore = 0) are left unchanged: the Glide config format is frozen.

Only the standard library is used.
"""

import argparse
import re
import sys

# (old section, old key) -> (new section, new key).
# Filled in one topic at a time, in the commits of the config reorganization.
MOVES = {
    # [calving]
    ('options', 'marine_margin'): ('calving', 'calving_law'),
    ('options', 'calving_init'): ('calving', 'calving_init'),
    ('options', 'calving_domain'): ('calving', 'calving_domain'),
    ('options', 'apply_calving_mask'): ('calving', 'apply_calving_mask'),
    ('options', 'remove_icebergs'): ('calving', 'remove_icebergs'),
    ('options', 'remove_isthmuses'): ('calving', 'remove_isthmuses'),
    ('options', 'limit_marine_cliffs'): ('calving', 'limit_marine_cliffs'),
    ('ho_options', 'which_ho_calving_front'): ('calving', 'subgrid_calving_front'),
    ('ho_options', 'which_ho_calvingmip_domain'): ('calving', 'calvingmip_domain'),
    ('parameters', 'marine_limit'): ('calving', 'marine_limit'),
    ('parameters', 'calving_fraction'): ('calving', 'calving_fraction'),
    ('parameters', 'calving_timescale'): ('calving', 'calving_timescale'),
    ('parameters', 'calving_minthck'): ('calving', 'calving_minthck'),
    ('parameters', 'dthck_dx_cf'): ('calving', 'dthck_dx_cf'),
    ('parameters', 'thck_effective_min'): ('calving', 'thck_effective_min'),
    ('parameters', 'eigenconstant'): ('calving', 'eigenconstant'),
    ('parameters', 'tau_eigenconstant1'): ('calving', 'tau_eigenconstant1'),
    ('parameters', 'tau_eigenconstant2'): ('calving', 'tau_eigenconstant2'),
    ('parameters', 'stress_threshold'): ('calving', 'stress_threshold'),
    ('parameters', 'lateral_rate_min'): ('calving', 'lateral_rate_min'),
    ('parameters', 'effec_stress_min'): ('calving', 'effec_stress_min'),
    ('parameters', 'length_scale'): ('calving', 'calving_length_scale'),
    ('parameters', 'damage_constant'): ('calving', 'damage_constant'),
    ('parameters', 'taumax_cliff'): ('calving', 'taumax_cliff'),
    ('parameters', 'cliff_timescale'): ('calving', 'cliff_timescale'),
    ('parameters', 'calving_front_x'): ('calving', 'calving_front_x'),
    ('parameters', 'calving_front_y'): ('calving', 'calving_front_y'),
    ('parameters', 'calving_front_radius'): ('calving', 'calving_front_radius'),
    ('parameters', 'cf_advance_retreat_amplitude'): ('calving', 'cf_advance_retreat_amplitude'),
    ('parameters', 'cf_advance_retreat_period'): ('calving', 'cf_advance_retreat_period'),
    ('parameters', 'f_ground_threshold'): ('calving', 'f_ground_threshold'),
    # [removal]
    ('ho_options', 'remove_ice_caps'): ('removal', 'remove_ice_caps'),
    ('ho_options', 'force_retreat'): ('removal', 'force_retreat'),
    # [basal_physics]
    ('options', 'slip_coeff'): ('basal_physics', 'btrc_method'),
    ('ho_options', 'which_ho_babc'): ('basal_physics', 'sliding_law'),
    ('ho_options', 'use_c_space_factor'): ('basal_physics', 'use_c_space_factor'),
    ('ho_options', 'which_ho_beta_limit'): ('basal_physics', 'beta_limit_method'),
    ('ho_options', 'which_ho_powerlaw_c'): ('basal_physics', 'powerlaw_c_method'),
    ('ho_options', 'which_ho_coulomb_c'): ('basal_physics', 'coulomb_c_method'),
    ('ho_options', 'elevation_based_coulomb_c'): ('basal_physics', 'elevation_based_coulomb_c'),
    ('ho_options', 'which_ho_effecpress'): ('basal_physics', 'effecpress_method'),
    ('parameters', 'basal_tract_const'): ('basal_physics', 'btrc_const'),
    ('parameters', 'beta_grounded_min'): ('basal_physics', 'beta_grounded_min'),
    ('parameters', 'ho_beta_const'): ('basal_physics', 'beta_const'),
    ('parameters', 'ho_beta_small'): ('basal_physics', 'beta_small'),
    ('parameters', 'ho_beta_large'): ('basal_physics', 'beta_large'),
    ('parameters', 'powerlaw_c_const'): ('basal_physics', 'powerlaw_c_const'),
    ('parameters', 'powerlaw_c_max'): ('basal_physics', 'powerlaw_c_max'),
    ('parameters', 'powerlaw_c_min'): ('basal_physics', 'powerlaw_c_min'),
    ('parameters', 'powerlaw_m'): ('basal_physics', 'powerlaw_m'),
    ('parameters', 'coulomb_c_const'): ('basal_physics', 'coulomb_c_const'),
    ('parameters', 'coulomb_c_max'): ('basal_physics', 'coulomb_c_max'),
    ('parameters', 'coulomb_c_min'): ('basal_physics', 'coulomb_c_min'),
    ('parameters', 'coulomb_c_hi'): ('basal_physics', 'coulomb_c_hi'),
    ('parameters', 'coulomb_c_lo'): ('basal_physics', 'coulomb_c_lo'),
    ('parameters', 'coulomb_c_bed_hi'): ('basal_physics', 'coulomb_c_bed_hi'),
    ('parameters', 'coulomb_c_bed_lo'): ('basal_physics', 'coulomb_c_bed_lo'),
    ('parameters', 'schoof_n'): ('basal_physics', 'schoof_n'),
    ('parameters', 'schoof_gamma'): ('basal_physics', 'schoof_gamma'),
    ('parameters', 'schoof_p'): ('basal_physics', 'schoof_p'),
    ('parameters', 'beta_powerlaw_umax'): ('basal_physics', 'beta_powerlaw_umax'),
    ('parameters', 'zoet_iversion_ut'): ('basal_physics', 'zoet_iverson_ut'),
    ('parameters', 'zoet_iversion_nmax'): ('basal_physics', 'zoet_iverson_nmax'),
    ('parameters', 'friction_powerlaw_k'): ('basal_physics', 'friction_powerlaw_k'),
    ('parameters', 'coulomb_bump_max_slope'): ('basal_physics', 'coulomb_bump_max_slope'),
    ('parameters', 'coulomb_bump_wavelength'): ('basal_physics', 'coulomb_bump_wavelength'),
    ('parameters', 'p_ocean_penetration'): ('basal_physics', 'p_ocean_penetration'),
    ('parameters', 'pseudo_plastic_q'): ('basal_physics', 'pseudo_plastic_q'),
    ('parameters', 'pseudo_plastic_u0'): ('basal_physics', 'pseudo_plastic_u0'),
    ('parameters', 'pseudo_plastic_phimin'): ('basal_physics', 'pseudo_plastic_phimin'),
    ('parameters', 'pseudo_plastic_phimax'): ('basal_physics', 'pseudo_plastic_phimax'),
    ('parameters', 'pseudo_plastic_bedmin'): ('basal_physics', 'pseudo_plastic_bedmin'),
    ('parameters', 'pseudo_plastic_bedmax'): ('basal_physics', 'pseudo_plastic_bedmax'),
    # [grounding_line]
    ('ho_options', 'which_ho_ground'): ('grounding_line', 'grounding_line_param'),
    ('ho_options', 'which_ho_fground_no_glp'): ('grounding_line', 'fground_no_glp'),
    ('ho_options', 'which_ho_ground_bmlt'): ('grounding_line', 'grounding_line_bmlt'),
    ('ho_options', 'which_ho_flotation_function'): ('grounding_line', 'flotation_function'),
    # [basal_hydro] (existing section): the switch moves in; one key renamed in place
    ('ho_options', 'which_ho_bwat'): ('basal_hydro', 'bwat_method'),
    ('basal_hydro', 'ho_flux_routing_scheme'): ('basal_hydro', 'flux_routing_method'),
    # [climate]
    ('options', 'smb_input'): ('climate', 'smb_input'),
    ('options', 'smb_input_function'): ('climate', 'smb_input_function'),
    ('options', 'artm_input_function'): ('climate', 'artm_input_function'),
    ('options', 'enable_smb_anomaly'): ('climate', 'enable_smb_anomaly'),
    ('options', 'enable_artm_anomaly'): ('climate', 'enable_artm_anomaly'),
    ('options', 'enable_snow_anomaly'): ('climate', 'enable_snow_anomaly'),
    ('options', 'enable_precip_anomaly'): ('climate', 'enable_precip_anomaly'),
    ('options', 'overwrite_acab'): ('climate', 'overwrite_acab'),
    ('options', 'enable_acab_dthck_dt_correction'): ('climate', 'enable_acab_dthck_dt_correction'),
    ('ho_options', 'block_inception'): ('climate', 'block_inception'),
    ('parameters', 't_lapse'): ('climate', 't_lapse'),
    ('parameters', 'degree_factor'): ('climate', 'degree_factor'),
    ('parameters', 'tmlt'): ('climate', 'tmlt'),
    ('parameters', 'snow_threshold_min'): ('climate', 'snow_threshold_min'),
    ('parameters', 'snow_threshold_max'): ('climate', 'snow_threshold_max'),
    ('parameters', 'smb_factor'): ('climate', 'smb_factor'),
    ('parameters', 'smb_anomaly_tstart'): ('climate', 'smb_anomaly_tstart'),
    ('parameters', 'smb_anomaly_timescale'): ('climate', 'smb_anomaly_timescale'),
    ('parameters', 'overwrite_acab_value'): ('climate', 'overwrite_acab_value'),
    ('parameters', 'overwrite_acab_minthck'): ('climate', 'overwrite_acab_minthck'),
    ('parameters', 'artm_anomaly_const'): ('climate', 'artm_anomaly_const'),
    ('parameters', 'artm_anomaly_tstart'): ('climate', 'artm_anomaly_tstart'),
    ('parameters', 'artm_anomaly_timescale'): ('climate', 'artm_anomaly_timescale'),
    # [basal_melt]
    ('options', 'bmlt_float'): ('basal_melt', 'bmlt_float_method'),
    ('options', 'bmlt_float_thermal_forcing_param'): ('basal_melt', 'bmlt_float_thermal_forcing_param'),
    ('options', 'bmlt_float_init'): ('basal_melt', 'bmlt_float_init'),
    ('ho_options', 'which_ho_deltaT_ocn'): ('basal_melt', 'deltaT_ocn_method'),
    ('options', 'enable_bmlt_anomaly'): ('basal_melt', 'enable_bmlt_anomaly'),
    ('options', 'basal_mass_balance'): ('basal_melt', 'basal_mass_balance'),
    ('parameters', 'bmlt_float_factor'): ('basal_melt', 'bmlt_float_factor'),
    ('parameters', 'bmlt_cavity_h0'): ('basal_melt', 'bmlt_cavity_h0'),
    ('parameters', 'bmlt_anomaly_tstart'): ('basal_melt', 'bmlt_anomaly_tstart'),
    ('parameters', 'bmlt_anomaly_timescale'): ('basal_melt', 'bmlt_anomaly_timescale'),
    ('parameters', 'bmlt_float_omega'): ('basal_melt', 'bmlt_float_omega'),
    ('parameters', 'bmlt_float_h0'): ('basal_melt', 'bmlt_float_h0'),
    ('parameters', 'bmlt_float_z0'): ('basal_melt', 'bmlt_float_z0'),
    ('parameters', 'bmlt_float_const'): ('basal_melt', 'bmlt_float_const'),
    ('parameters', 'bmlt_float_xlim'): ('basal_melt', 'bmlt_float_xlim'),
    ('parameters', 'bmlt_float_depth_meltmax'): ('basal_melt', 'bmlt_float_depth_meltmax'),
    ('parameters', 'bmlt_float_depth_frzmax'): ('basal_melt', 'bmlt_float_depth_frzmax'),
    ('parameters', 'bmlt_float_depth_zmeltmax'): ('basal_melt', 'bmlt_float_depth_zmeltmax'),
    ('parameters', 'bmlt_float_depth_zmelt0'): ('basal_melt', 'bmlt_float_depth_zmelt0'),
    ('parameters', 'bmlt_float_depth_zfrzmax'): ('basal_melt', 'bmlt_float_depth_zfrzmax'),
    ('parameters', 'bmlt_float_depth_meltmin'): ('basal_melt', 'bmlt_float_depth_meltmin'),
    ('parameters', 'bmlt_float_depth_zmeltmin'): ('basal_melt', 'bmlt_float_depth_zmeltmin'),
    # [ocean_data]
    ('options', 'ocean_data_domain'): ('ocean_data', 'ocean_data_domain'),
    ('options', 'ocean_data_extrapolate'): ('ocean_data', 'ocean_data_extrapolate'),
    ('parameters', 'gamma0'): ('ocean_data', 'gamma0'),
    ('parameters', 'thermal_forcing_basin_min'): ('ocean_data', 'thermal_forcing_basin_min'),
    ('parameters', 'thermal_forcing_basin_max'): ('ocean_data', 'thermal_forcing_basin_max'),
    ('parameters', 'thermal_forcing_anomaly'): ('ocean_data', 'thermal_forcing_anomaly'),
    ('parameters', 'thermal_forcing_anomaly_tstart'): ('ocean_data', 'thermal_forcing_anomaly_tstart'),
    ('parameters', 'thermal_forcing_anomaly_timescale'): ('ocean_data', 'thermal_forcing_anomaly_timescale'),
    ('parameters', 'thermal_forcing_anomaly_basin'): ('ocean_data', 'thermal_forcing_anomaly_basin'),
    # [velocity]
    ('ho_options', 'which_ho_approx'): ('velocity', 'stokes_approx'),
    ('ho_options', 'diva_slope_correction'): ('velocity', 'diva_slope_correction'),
    ('ho_options', 'which_ho_efvs'): ('velocity', 'efvs_method'),
    ('ho_options', 'which_ho_gradient'): ('velocity', 'surface_gradient_method'),
    ('ho_options', 'which_ho_gradient_margin'): ('velocity', 'surface_gradient_margin_method'),
    ('parameters', 'efvs_constant'): ('velocity', 'efvs_constant'),
    ('parameters', 'effstrain_min'): ('velocity', 'effstrain_min'),
    ('parameters', 'max_slope'): ('velocity', 'max_slope'),
    # [velocity_solver]
    ('ho_options', 'which_ho_nonlinear'): ('velocity_solver', 'nonlinear_method'),
    ('ho_options', 'which_ho_resid'): ('velocity_solver', 'nonlinear_residual'),
    ('ho_options', 'glissade_maxiter'): ('velocity_solver', 'nonlinear_maxiters'),
    ('ho_options', 'which_ho_sparse'): ('velocity_solver', 'linear_solver'),
    ('ho_options', 'which_ho_precond'): ('velocity_solver', 'preconditioner'),
    ('ho_options', 'linear_solve_ncheck'): ('velocity_solver', 'linear_solve_ncheck'),
    ('ho_options', 'linear_maxiters'): ('velocity_solver', 'linear_maxiters'),
    ('ho_options', 'linear_tolerance'): ('velocity_solver', 'linear_tolerance'),
    ('ho_options', 'which_ho_assemble_beta'): ('velocity_solver', 'assemble_beta'),
    ('ho_options', 'which_ho_assemble_taud'): ('velocity_solver', 'assemble_taud'),
    ('ho_options', 'which_ho_assemble_bfric'): ('velocity_solver', 'assemble_bfric'),
    ('ho_options', 'which_ho_assemble_lateral'): ('velocity_solver', 'assemble_lateral'),
    # [advection]
    ('options', 'evolution'): ('advection', 'advection_method'),
    ('ho_options', 'which_ho_vertical_remap'): ('advection', 'vertical_remap_order'),
    ('ho_options', 'which_ho_ice_age'): ('advection', 'compute_ice_age'),
    # [flow_law]
    ('options', 'flow_law'): ('flow_law', 'flow_law'),
    ('options', 'damage_flwa_feedback'): ('flow_law', 'damage_flwa_feedback'),
    ('parameters', 'n_glen'): ('flow_law', 'n_glen'),
    ('parameters', 'flow_factor_ground'): ('flow_law', 'flow_factor_ground'),
    ('parameters', 'flow_factor_float'): ('flow_law', 'flow_factor_float'),
    ('parameters', 'default_flwa'): ('flow_law', 'default_flwa'),
    # [temperature]
    ('options', 'temperature'): ('temperature', 'temperature_method'),
    ('options', 'temp_init'): ('temperature', 'temp_init'),
    ('ho_options', 'which_ho_disp'): ('temperature', 'dissipation_method'),
    ('ho_options', 'which_ho_thermal_timestep'): ('temperature', 'thermal_timestep'),
    ('parameters', 'pmp_offset'): ('temperature', 'pmp_offset'),
    ('parameters', 'pmp_threshold'): ('temperature', 'pmp_threshold'),
    # [inversion]
    ('parameters', 'inversion_thck_flotation_buffer'): ('inversion', 'inversion_thck_flotation_buffer'),
    ('parameters', 'inversion_thck_threshold'): ('inversion', 'inversion_thck_threshold'),
    ('parameters', 'inversion_damping_factor'): ('inversion', 'inversion_damping_factor'),
    ('parameters', 'inversion_babc_timescale'): ('inversion', 'inversion_babc_timescale'),
    ('parameters', 'inversion_babc_thck_scale'): ('inversion', 'inversion_babc_thck_scale'),
    ('parameters', 'inversion_babc_length_scale'): ('inversion', 'inversion_babc_length_scale'),
    ('parameters', 'inversion_babc_relax_factor'): ('inversion', 'inversion_babc_relax_factor'),
    ('parameters', 'inversion_deltaT_ocn_timescale'): ('inversion', 'inversion_deltaT_ocn_timescale'),
    ('parameters', 'inversion_deltaT_ocn_thck_scale'): ('inversion', 'inversion_deltaT_ocn_thck_scale'),
    ('parameters', 'inversion_deltaT_ocn_length_scale'): ('inversion', 'inversion_deltaT_ocn_length_scale'),
    ('parameters', 'inversion_deltaT_ocn_temp_scale'): ('inversion', 'inversion_deltaT_ocn_temp_scale'),
    ('parameters', 'inversion_deltaT_basin_relax'): ('inversion', 'inversion_deltaT_basin_relax'),
    ('parameters', 'inversion_basin_flotation_threshold'): ('inversion', 'inversion_basin_flotation_threshold'),
    ('parameters', 'inversion_basin_mass_correction'): ('inversion', 'inversion_basin_mass_correction'),
    ('parameters', 'inversion_basin_number_mass_correction'): ('inversion', 'inversion_basin_number_mass_correction'),
    # [geothermal]
    ('options', 'gthf'): ('geothermal', 'geothermal_method'),
    ('parameters', 'geothermal'): ('geothermal', 'geothermal_flux'),
    # [options]
    ('parameters', 'log_level'): ('options', 'log_level'),
    # [grid]
    ('options', 'read_lat_lon'): ('grid', 'read_lat_lon'),
    ('options', 'sigma'): ('grid', 'sigma_method'),
    # [numerics]
    ('ho_options', 'compute_blocks'): ('numerics', 'compute_blocks'),
    ('ho_options', 'reproducible_sums'): ('numerics', 'reproducible_sums'),
    ('parameters', 'ice_limit'): ('numerics', 'ice_limit'),
    ('parameters', 'ice_limit_temp'): ('numerics', 'ice_limit_temp'),
    ('parameters', 'thck_gradient_ramp'): ('numerics', 'thck_gradient_ramp'),
    ('parameters', 'periodic_offset_ew'): ('numerics', 'periodic_offset_ew'),
    ('parameters', 'periodic_offset_ns'): ('numerics', 'periodic_offset_ns'),
    # [constants]
    ('parameters', 'rhoi'): ('constants', 'rhoi'),
    ('parameters', 'rhoo'): ('constants', 'rhoo'),
    ('parameters', 'grav'): ('constants', 'grav'),
    ('parameters', 'shci'): ('constants', 'shci'),
    ('parameters', 'lhci'): ('constants', 'lhci'),
    ('parameters', 'trpt'): ('constants', 'trpt'),
    # [initialization]
    ('options', 'adjust_input_thickness'): ('initialization', 'adjust_input_thickness'),
    ('options', 'smooth_input_topography'): ('initialization', 'smooth_input_topography'),
    ('options', 'smooth_input_usrf'): ('initialization', 'smooth_input_usrf'),
    ('options', 'adjust_input_topography'): ('initialization', 'adjust_input_topography'),
    ('parameters', 'adjust_topg_xmin'): ('initialization', 'adjust_topg_xmin'),
    ('parameters', 'adjust_topg_xmax'): ('initialization', 'adjust_topg_xmax'),
    ('parameters', 'adjust_topg_ymin'): ('initialization', 'adjust_topg_ymin'),
    ('parameters', 'adjust_topg_ymax'): ('initialization', 'adjust_topg_ymax'),
    ('parameters', 'adjust_topg_no_adjust'): ('initialization', 'adjust_topg_no_adjust'),
    ('parameters', 'adjust_topg_max_adjust'): ('initialization', 'adjust_topg_max_adjust'),
    ('parameters', 'adjust_topg_delta'): ('initialization', 'adjust_topg_delta'),
    # [diagnostics]
    ('time', 'dt_diag'): ('diagnostics', 'dt_diag'),
    ('time', 'ndiag'): ('diagnostics', 'ndiag'),
    ('time', 'idiag'): ('diagnostics', 'idiag'),
    ('time', 'jdiag'): ('diagnostics', 'jdiag'),
    ('options', 'dm_dt_diag'): ('diagnostics', 'dm_dt_diag'),
    ('options', 'diag_minthck'): ('diagnostics', 'diag_minthck'),
    # [isostasy], [lateral_melt], [glaciers]: section switches move in; keys renamed in place
    ('options', 'isostasy'): ('isostasy', 'isostasy'),
    ('isostasy', 'lithosphere'): ('isostasy', 'lithosphere_method'),
    ('isostasy', 'asthenosphere'): ('isostasy', 'asthenosphere_method'),
    ('isostasy', 'which_relaxed'): ('isostasy', 'input_topg_method'),
    ('options', 'which_lateral_melt'): ('lateral_melt', 'lateral_melt_method'),
    ('ho_options', 'enable_glaciers'): ('glaciers', 'enable_glaciers'),
    ('glaciers', 'set_mu_star'): ('glaciers', 'mu_star_method'),
    ('glaciers', 'set_alpha_snow'): ('glaciers', 'alpha_snow_method'),
    ('glaciers', 'set_powerlaw_c'): ('glaciers', 'glacier_powerlaw_c_method'),
    ('glaciers', 'snow_calc'): ('glaciers', 'snow_method'),
    ('glaciers', 'tmlt'): ('glaciers', 'glacier_tmlt'),
    ('glaciers', 'snow_threshold_min'): ('glaciers', 'glacier_snow_threshold_min'),
    ('glaciers', 'snow_threshold_max'): ('glaciers', 'glacier_snow_threshold_max'),
}

# (section, key) -> note: keys that Glissade no longer reads.  The converter deletes them and reports the
# note, which names the replacement if there is one.
_FLOW_ENH = 'removed with the flow-enhancement-factor inversion'
_PLUME = 'removed with the old plume model'
REMOVED = {
    ('ho_options', 'which_ho_flow_enhancement_factor'): 'removed (only the default behavior remains)',
    ('parameters', 'inversion_flow_enhancement_timescale'): _FLOW_ENH,
    ('parameters', 'inversion_flow_enhancement_velo_scale'): _FLOW_ENH,
    ('parameters', 'inversion_flow_enhancement_thck_scale'): _FLOW_ENH,
    ('parameters', 'inversion_flow_enhancement_length_scale'): _FLOW_ENH,
    ('parameters', 'inversion_flow_enhancement_relax_factor'): _FLOW_ENH,
    ('parameters', 'flwa_basal'): 'removed (unused; [basal_hydro] flwa_basal is a different key)',
    ('parameters', 'T0'): _PLUME, ('parameters', 'Tbot'): _PLUME,
    ('parameters', 'S0'): _PLUME, ('parameters', 'Sbot'): _PLUME,
    ('parameters', 'gammaT'): _PLUME, ('parameters', 'gammaS'): _PLUME,
    ('ho_options', 'deltaT_ocn_extrapolate'): 'removed (deltaT_ocn is no longer extrapolated)',
    ('isostasy', 'lithosphere_period'): 'removed; see load_update_interval in [isostasy]',
    ('options', 'enable_acab_anomaly'): 'removed; see enable_smb_anomaly in [climate]',
}

# Order of the Glissade sections in a converted config file: the order in which they are written
# to the log (glissade_printconfig), with [sigma] after [grid].  Sections not listed here ([CF ...]
# and any others) follow, in their original order.
SECTION_ORDER = [
    'grid', 'sigma', 'time', 'diagnostics', 'options', 'numerics', 'constants', 'initialization',
    'velocity', 'velocity_solver', 'advection', 'temperature', 'flow_law',
    'basal_physics', 'basal_hydro', 'grounding_line', 'basal_melt', 'calving', 'lateral_melt', 'removal',
    'grid_atm', 'climate', 'grid_ocn', 'ocean_data', 'geothermal', 'isostasy', 'inversion', 'glaciers',
]

# CESM user_nl_cism: old key -> new key.  Kept separate from MOVES because a few key names occur in
# more than one CISM section, while the CESM namelist has one entry per name.
USER_NL_RENAMES = {
    # [calving]
    'marine_margin': 'calving_law',
    'which_ho_calving_front': 'subgrid_calving_front',
    'which_ho_calvingmip_domain': 'calvingmip_domain',
    'length_scale': 'calving_length_scale',
    # [basal_physics]
    'which_ho_babc': 'sliding_law',
    'which_ho_beta_limit': 'beta_limit_method',
    'which_ho_powerlaw_c': 'powerlaw_c_method',
    'which_ho_coulomb_c': 'coulomb_c_method',
    'which_ho_effecpress': 'effecpress_method',
    'basal_tract_const': 'btrc_const',
    'ho_beta_const': 'beta_const',
    # [grounding_line]
    'which_ho_ground': 'grounding_line_param',
    'which_ho_ground_bmlt': 'grounding_line_bmlt',
    'which_ho_flotation_function': 'flotation_function',
    # [basal_hydro]
    'which_ho_bwat': 'bwat_method',
    'ho_flux_routing_scheme': 'flux_routing_method',
    # [basal_melt]
    'bmlt_float': 'bmlt_float_method',
    'which_ho_deltaT_ocn': 'deltaT_ocn_method',
    # [velocity], [velocity_solver]
    'which_ho_approx': 'stokes_approx',
    'which_ho_efvs': 'efvs_method',
    'which_ho_gradient': 'surface_gradient_method',
    'which_ho_gradient_margin': 'surface_gradient_margin_method',
    'which_ho_nonlinear': 'nonlinear_method',
    'which_ho_resid': 'nonlinear_residual',
    'glissade_maxiter': 'nonlinear_maxiters',
    'which_ho_sparse': 'linear_solver',
    'which_ho_precond': 'preconditioner',
    'which_ho_assemble_beta': 'assemble_beta',
    'which_ho_assemble_taud': 'assemble_taud',
    'which_ho_assemble_bfric': 'assemble_bfric',
    'which_ho_assemble_lateral': 'assemble_lateral',
    # [advection], [temperature]
    'evolution': 'advection_method',
    'which_ho_vertical_remap': 'vertical_remap_order',
    'which_ho_ice_age': 'compute_ice_age',
    'temperature': 'temperature_method',
    'which_ho_disp': 'dissipation_method',
    'which_ho_thermal_timestep': 'thermal_timestep',
    # [geothermal]
    'gthf': 'geothermal_method',
    'geothermal': 'geothermal_flux',
    # [grid]
    'sigma': 'sigma_method',
    # [isostasy], [lateral_melt], [glaciers] (in CESM, tmlt and snow_threshold_* are the glacier keys)
    'which_lateral_melt': 'lateral_melt_method',
    'which_relaxed': 'input_topg_method',
    'lithosphere': 'lithosphere_method',
    'asthenosphere': 'asthenosphere_method',
    'set_mu_star': 'mu_star_method',
    'set_alpha_snow': 'alpha_snow_method',
    'set_powerlaw_c': 'glacier_powerlaw_c_method',
    'snow_calc': 'snow_method',
    'tmlt': 'glacier_tmlt',
    'snow_threshold_min': 'glacier_snow_threshold_min',
    'snow_threshold_max': 'glacier_snow_threshold_max',
    # [inversion]: wrapper names that CISM never read (the wrapper now uses the CISM names)
    'deltaT_basin_relax': 'inversion_deltaT_basin_relax',
    'inversion_bmlt_basin_flotation_threshold': 'inversion_basin_flotation_threshold',
    'inversion_bmlt_basin_mass_correction': 'inversion_basin_mass_correction',
    'inversion_bmlt_basin_number_mass_correction': 'inversion_basin_number_mass_correction',
}

# Match the CISM reader (glimmer_config.F90): a section header is a line starting with '[' in
# column 1, named by the text before the first ']' (anything after it is ignored); a key line
# starts in column 1 (lines starting with a space are ignored by CISM, so they are never moved).
SECTION_RE = re.compile(r'^\[([^\]]*)\]')
KEY_RE = re.compile(r'^()([A-Za-z_][A-Za-z0-9_]*)(\s*[=:].*)$')
COMMENT_RE = re.compile(r'^\s*[#;!]')
NAMELIST_KEY_RE = re.compile(r'^(\s*)([A-Za-z_][A-Za-z0-9_]*)(\s*=.*)$')


class Section:
    def __init__(self, name, header):
        self.name = name        # None for lines before the first section
        self.header = header    # the original header line
        self.lines = []


def parse(text):
    """Split a config file into sections, keeping every line."""
    sections = [Section(None, None)]
    for line in text.split('\n'):
        m = SECTION_RE.match(line)
        if m:
            sections.append(Section(m.group(1).strip(), line))
        else:
            sections[-1].lines.append(line)
    return sections


def render(sections):
    out = []
    for s in sections:
        if s.header is not None:
            out.append(s.header)
        out.extend(s.lines)
    return '\n'.join(out)


def has_keys(section):
    return any(KEY_RE.match(l) and not COMMENT_RE.match(l) for l in section.lines)


def is_glide(sections):
    """True if [options] sets dycore = 0 (Glide)."""
    for s in sections:
        if s.name == 'options':
            for line in s.lines:
                m = KEY_RE.match(line)
                if m and not COMMENT_RE.match(line) and m.group(2) == 'dycore':
                    value = re.split(r'[#;!]', m.group(3).lstrip(' \t=:'))[0].strip()
                    return value == '0'
    return False


def reorder_sections(sections):
    """Put the sections in SECTION_ORDER, followed by the other sections in their original order.

    Returns True if the order changed.  If it did not, nothing is changed (not even whitespace).
    Comment lines directly above a section header, separated from the previous section's content
    by a blank line, move with the section below them (e.g. '# I/O configuration follows' above
    [CF default]).  Lines before the first section stay at the top of the file.
    """
    head, body = sections[:1], sections[1:]          # sections[0] holds the lines before the first header
    rank = {name: i for i, name in enumerate(SECTION_ORDER)}
    order = sorted(range(len(body)), key=lambda i: (rank.get(body[i].name, len(SECTION_ORDER)), i))
    if order == list(range(len(body))):
        return False
    # detach the comment block above each header from the end of the preceding section
    lead = [[] for _ in body]
    for i in range(1, len(body)):
        prev = body[i - 1].lines
        n = len(prev)
        while n > 0 and prev[n - 1].strip() and COMMENT_RE.match(prev[n - 1]):
            n -= 1
        if n < len(prev) and (n == 0 or not prev[n - 1].strip()):
            lead[i], body[i - 1].lines = prev[n:], prev[:n]
    new = []
    for k, i in enumerate(order):
        s = body[i]
        if k < len(order) - 1 and (not s.lines or s.lines[-1].strip()):
            s.lines.append('')                        # keep a blank line between sections
        new.append((lead[i], s))
    sections[1:] = []
    for comments, s in new:
        if comments:
            # the comment block goes above the header: store it at the end of the section before it
            target = sections[-1]
            target.lines.extend(comments)
        sections.append(s)
    # a file that ended with a newline still does
    last = sections[-1].lines
    if body and order[-1] != len(body) - 1 and (not last or last[-1] != ''):
        last.append('')
    return True


def convert_config(text, moves=None, removed=None):
    """Return (new_text, list of changes) for a CISM config file.  Glide configs are returned unchanged."""
    moves = MOVES if moves is None else moves
    removed = REMOVED if removed is None else removed
    sections = parse(text)
    if is_glide(sections):
        return text, []
    changes = []
    moved = {}   # new section -> list of lines to add
    emptied = set()   # sections that lost keys

    for s in sections:
        if s.name is None:
            continue
        kept = []
        for line in s.lines:
            m = KEY_RE.match(line)
            if m and not COMMENT_RE.match(line) and (s.name, m.group(2)) in moves:
                new_sec, new_key = moves[(s.name, m.group(2))]
                # take the comment lines directly above this key along with it
                carried = []
                while kept and kept[-1].strip() and COMMENT_RE.match(kept[-1]):
                    carried.insert(0, kept.pop())
                moved.setdefault(new_sec, []).extend(carried + [new_key + m.group(3)])
                emptied.add(s.name)
                changes.append(f'[{s.name}] {m.group(2)} -> [{new_sec}] {new_key}')
            elif m and not COMMENT_RE.match(line) and (s.name, m.group(2)) in removed:
                emptied.add(s.name)
                changes.append(f'deleted [{s.name}] {m.group(2)}: {removed[(s.name, m.group(2))]}')
            else:
                kept.append(line)
        s.lines = kept

    for new_sec, lines in moved.items():
        target = next((s for s in sections if s.name == new_sec), None)
        if target is None:
            target = Section(new_sec, f'[{new_sec}]')
            target.lines = ['']
            # insert before the first [CF ...] section, keeping SECTION_ORDER among new sections
            pos = next((i for i, s in enumerate(sections) if s.name and s.name.startswith('CF ')), len(sections))
            sections.insert(pos, target)
            changes.append(f'new section [{new_sec}]')
        existing = {KEY_RE.match(l).group(2) for l in target.lines if KEY_RE.match(l) and not COMMENT_RE.match(l)}
        for l in lines:
            m = KEY_RE.match(l)
            if m and not COMMENT_RE.match(l) and m.group(2) in existing:
                raise ValueError(f'[{new_sec}] already has a key {m.group(2)}')
        # add before the trailing blank lines of the target section
        n = len(target.lines)
        while n > 0 and not target.lines[n - 1].strip():
            n -= 1
        target.lines[n:n] = lines

    # remove old sections that have no keys left
    for s in list(sections):
        if s.name in emptied and s.name not in moved and not has_keys(s):
            sections.remove(s)
            changes.append(f'removed empty section [{s.name}]')

    # put the sections in the standard order
    if reorder_sections(sections):
        changes.append('sections reordered')

    return render(sections), changes


def convert_user_nl(text, renames=None):
    """Return (new_text, list of changes) for a CESM user_nl_cism file."""
    renames = USER_NL_RENAMES if renames is None else renames
    out, changes = [], []
    for line in text.split('\n'):
        m = NAMELIST_KEY_RE.match(line)
        if m and not COMMENT_RE.match(line) and m.group(2) in renames:
            out.append(m.group(1) + renames[m.group(2)] + m.group(3))
            changes.append(f'{m.group(2)} -> {renames[m.group(2)]}')
        else:
            out.append(line)
    return '\n'.join(out), changes


def is_config(text):
    return any(SECTION_RE.match(l) for l in text.split('\n'))


def self_test():
    moves = {('options', 'marine_margin'): ('calving', 'calving_law'),
             ('parameters', 'calving_timescale'): ('calving', 'calving_timescale'),
             ('ho_options', 'which_ho_babc'): ('basal_physics', 'sliding_law')}
    old = '\n'.join([
        '[grid]', 'ewn = 31', '',
        '[options]', 'dycore = 2', '# calving: 4 = thickness threshold', 'marine_margin = 4   # see docs', '',
        '[ho_options]', 'which_ho_babc = 9', '',
        '[parameters]', 'calving_timescale = 1.0', 'ice_limit = 1.', '',
        '[CF output]', 'name = out.nc', ''])
    new, changes = convert_config(old, moves)
    assert '[calving]\n# calving: 4 = thickness threshold\ncalving_law = 4   # see docs\ncalving_timescale = 1.0\n\n' in new, new
    assert '[basal_physics]\nsliding_law = 9\n\n' in new, new
    assert '[ho_options]' not in new, new                      # emptied section removed
    assert new.index('[basal_physics]') < new.index('[calving]') < new.index('[CF output]'), new
    assert 'dycore = 2' in new and 'ice_limit = 1.' in new
    again, changes2 = convert_config(new, moves)                # idempotent
    assert again == new and not changes2, changes2
    nl, c = convert_user_nl("which_ho_babc = 9\n! comment\nice_limit = 1.\n", {'which_ho_babc': 'sliding_law'})
    assert nl == "sliding_law = 9\n! comment\nice_limit = 1.\n" and c == ['which_ho_babc -> sliding_law'], nl
    assert convert_config(old, {})[0] == old                    # empty table: no change
    commented = old.replace('[ho_options]', '[ho_options]   # Higher order options').replace(
        'dycore = 2', 'dycore = 2\n marine_margin = 1')
    new2, _ = convert_config(commented, moves)                  # header with a trailing comment
    assert '[basal_physics]\nsliding_law = 9\n' in new2 and '[ho_options]' not in new2, new2
    assert ' marine_margin = 1' in new2                         # indented key (ignored by CISM) left alone
    glide = old.replace('dycore = 2', 'dycore = 0   # Glide')
    assert convert_config(glide, moves) == (glide, [])          # Glide config: no change

    # section order: comment blocks above a header (after a blank line) move with the section below
    unordered = '\n'.join([
        '# dome test', '[options]', 'dycore = 2', '',
        '# basal friction', '[basal_physics]', 'sliding_law = 4', '#beta_large = 1.e10',
        '[grid]', 'ewn = 31', '',
        '####', '# I/O configuration follows', '####', '[CF default]', 'title = test', '',
        '[time]', 'dt = 1.', ''])
    ordered, ch = convert_config(unordered, {})
    assert ch == ['sections reordered'], ch
    assert ordered == '\n'.join([
        '# dome test', '[grid]', 'ewn = 31', '',
        '[time]', 'dt = 1.', '',
        '[options]', 'dycore = 2', '',
        '# basal friction', '[basal_physics]', 'sliding_law = 4', '#beta_large = 1.e10', '',
        '####', '# I/O configuration follows', '####', '[CF default]', 'title = test', '']), ordered
    assert convert_config(ordered, {}) == (ordered, [])          # already in order: unchanged

    # removed keys are deleted and reported; a section left with no keys is removed
    stale = '\n'.join(['[grid]', 'ewn = 31', '', '[ho_options]', 'deltaT_ocn_extrapolate = .true.', '',
                        '[isostasy]', 'isostasy = 1', 'lithosphere_period = 10.   # yr', ''])
    clean, ch = convert_config(stale, {})
    assert clean == '\n'.join(['[grid]', 'ewn = 31', '', '[isostasy]', 'isostasy = 1', '']), clean
    assert ch == ['deleted [ho_options] deltaT_ocn_extrapolate: removed (deltaT_ocn is no longer extrapolated)',
                  'deleted [isostasy] lithosphere_period: removed; see load_update_interval in [isostasy]',
                  'removed empty section [ho_options]'], ch
    print('self-test passed')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('files', nargs='*')
    ap.add_argument('-o', '--output')
    ap.add_argument('--in-place', action='store_true')
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--user-nl', action='store_true', help='treat the input as a CESM user_nl_cism file')
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    if a.self_test:
        self_test()
        return 0
    if not a.files:
        ap.error('no input files')
    if a.output and len(a.files) != 1:
        ap.error('-o needs exactly one input file')
    any_change = False
    for f in a.files:
        text = open(f).read()
        if a.user_nl or not is_config(text):
            new, changes = convert_user_nl(text)
        else:
            new, changes = convert_config(text)
        any_change |= bool(changes)
        glide = not (a.user_nl or not is_config(text)) and is_glide(parse(text))
        if a.check:
            print(f'{f}: ' + ('Glide config, left unchanged' if glide else '; '.join(changes) if changes else 'no changes'))
            continue
        dest = a.output or (f if a.in_place else f + '.converted')
        if changes or dest != f:
            open(dest, 'w').write(new)
        print(f'{f}: {len(changes)} changes' + (f', written to {dest}' if dest != f else ''))
    return 1 if (a.check and any_change) else 0


if __name__ == '__main__':
    sys.exit(main())
