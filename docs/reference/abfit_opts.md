# Return a list of options for an ABfit analysis.

Return a list of options for an ABfit analysis.

## Usage

``` r
abfit_opts(
  init_damping = 5,
  maxiters = 1024,
  max_shift_pre = 0.078,
  max_shift_fine = NULL,
  max_damping = 15,
  max_phase = 360,
  lambda = NULL,
  ppm_left = 4,
  ppm_right = 0.2,
  zp = TRUE,
  bl_ed_pppm = 2,
  auto_bl_flex = TRUE,
  bl_comps_pppm = 15,
  adaptive_bl_comps_pppm = FALSE,
  export_sp_fit = FALSE,
  max_asym = 0.25,
  max_basis_shift = 0.0078,
  max_basis_damping = 2,
  maxiters_pre = 1000,
  algo_pre = "NLOPT_LN_NELDERMEAD",
  min_bl_ed_pppm = NULL,
  max_bl_ed_pppm = 7,
  auto_bl_flex_n = 20,
  pre_fit_bl_ed_pppm = 1,
  remove_lip_mm_prefit = FALSE,
  pre_align = TRUE,
  max_pre_align_shift = 0.1,
  pre_align_ref_freqs = c(2.01, 3.03, 3.22),
  noise_region = c(-0.5, -2.5),
  optimal_smooth_criterion = "maic",
  aic_smoothing_factor = 5,
  anal_jac = TRUE,
  pre_fit_ppm_left = 4,
  pre_fit_ppm_right = 1.8,
  phi1_optim = FALSE,
  phi1_init = 0,
  max_dphi1 = 0.2,
  max_basis_shift_broad = NULL,
  max_basis_damping_broad = NULL,
  ahat_calc_method = "lh_pnnls",
  prefit_phase_search = TRUE,
  freq_reg = NULL,
  freq_reg_naa = NULL,
  lb_reg = NULL,
  asym_reg = NULL,
  output_all_paras = FALSE,
  output_all_paras_raw = FALSE,
  input_paras_raw = NULL,
  optim_lw_only = FALSE,
  optim_lw_only_limit = 20,
  lb_init = 0.001,
  lb_init_approx_fit = FALSE,
  zf_offset = NULL,
  broad_asym = TRUE,
  lb_reg_broad = NULL,
  broad_glb = TRUE
)
```

## Arguments

- init_damping:

  initial value of the Gaussian global damping parameter (Hz). Very
  poorly shimmed or high field data may benefit from a larger value.

- maxiters:

  The maximum number of iterations to run for the detailed fit.

- max_shift_pre:

  The maximum allowable global shift to be applied in the approximate
  (pre-fit) phases of analysis (ppm).

- max_shift_fine:

  The maximum allowable global shift to be applied in the detailed fit
  phase of analysis (ppm).

- max_damping:

  maximum permitted value of the global damping parameter (Hz).

- max_phase:

  the maximum absolute permitted value of the global zero-order phase
  term (degrees). Note, the prefit_phase_search option is not
  constrained by this term.

- lambda:

  manually set the the baseline smoothness parameter.

- ppm_left:

  downfield frequency limit for the fitting range (ppm).

- ppm_right:

  upfield frequency limit for the fitting range (ppm).

- zp:

  zero pad the data to twice the original length before fitting.

- bl_ed_pppm:

  manually set the the baseline smoothness parameter (ED per ppm).

- auto_bl_flex:

  automatically determine the level of baseline smoothness.

- bl_comps_pppm:

  spline basis density (signals per ppm).

- adaptive_bl_comps_pppm:

  adjust the spline basis density in the detailed fit phase, based on
  the required level of smoothness, to reduce computation time.

- export_sp_fit:

  add the fitted spline functions to the fit result.

- max_asym:

  maximum allowable value of the asymmetry parameter.

- max_basis_shift:

  maximum allowable frequency shift for individual basis signals (ppm).

- max_basis_damping:

  maximum allowable Lorentzian damping factor for individual basis
  signals (Hz).

- maxiters_pre:

  maximum iterations for the coarse (pre-)fit.

- algo_pre:

  optimisation method for the coarse (pre-)fit.

- min_bl_ed_pppm:

  minimum value for the candidate baseline flexibility analyses (ED per
  ppm).

- max_bl_ed_pppm:

  minimum value for the candidate baseline flexibility analyses (ED per
  ppm).

- auto_bl_flex_n:

  number of candidate baseline analyses to perform.

- pre_fit_bl_ed_pppm:

  level of baseline flexibility to use in the coarse fitting stage of
  the algorithm (ED per ppm).

- remove_lip_mm_prefit:

  remove broad signals in the coarse fitting stage of the algorithm.

- pre_align:

  perform a pre-alignment step before coarse fitting.

- max_pre_align_shift:

  maximum allowable shift in the pre-alignment step (ppm).

- pre_align_ref_freqs:

  a vector of prominent spectral frequencies used in the pre-alignment
  step (ppm).

- noise_region:

  spectral region to estimate the noise level (ppm).

- optimal_smooth_criterion:

  method to determine the optimal smoothness.

- aic_smoothing_factor:

  modification factor for the AIC calculation. Larger values result in
  less flexible baselines.

- anal_jac:

  use a analytical approximation to the jacobian in the detailed fitting
  stage.

- pre_fit_ppm_left:

  downfield frequency limit for the fitting range in the coarse fitting
  stage of the algorithm (ppm).

- pre_fit_ppm_right:

  upfield frequency limit for the fitting range in the coarse fitting
  stage of the algorithm (ppm).

- phi1_optim:

  apply and optimise a frequency dependant phase term.

- phi1_init:

  initial value for the frequency dependant phase term (ms).

- max_dphi1:

  maximum allowable change from the initial frequency dependant phase
  term (ms).

- max_basis_shift_broad:

  maximum allowable shift for broad signals in the basis (ppm).
  Determined based on their name beginning with Lip or MM. The default
  value is set to max_basis_shift.

- max_basis_damping_broad:

  maximum allowable Lorentzian damping for broad signals in the basis
  (Hz). Determined based on their name beginning with Lip or MM. The
  default value is set to max_basis_damping.

- ahat_calc_method:

  method to calculate the metabolite amplitudes. May be one of:
  "lh_pnnls" or "ls".

- prefit_phase_search:

  perform a 1D search for the optimal phase in the prefit stage of the
  algorithm.

- freq_reg:

  frequency shift parameter.

- freq_reg_naa:

  frequency shift parameter for NAA and NAAG.

- lb_reg:

  individual line broadening parameter.

- asym_reg:

  lineshape asymmetry parameter.

- output_all_paras:

  include more fitting parameters in the fit table, e.g. individual
  shift and damping factors for each basis set element.

- output_all_paras_raw:

  include raw fitting parameters in the fit table. For advanced
  diagnostic use only.

- input_paras_raw:

  input raw fitting parameters. For advanced diagnostic use only.

- optim_lw_only:

  optimize the global line-broadening term only.

- optim_lw_only_limit:

  limits for the line-breading term as a percentage of the starting
  value when optim_lw_only is TRUE.

- lb_init:

  initial Lorentzian line broadening value (in Hz) for the individual
  basis signals. Setting to 0 will clash with the minimum allowable
  value (eg hard constraint) during the detailed fit.

- lb_init_approx_fit:

  apply lb_init to the basis during the approximate iterative fit.

- zf_offset:

  offset in number of data points from the end of the FID to zero-fill.
  Default is NULL and will automatically set this to 50 points when the
  FID distortion flag is set for the mrs_data.

- broad_asym:

  apply asymmetric lineshape parameter to broad signals with names
  starting with Lip or MM.

- lb_reg_broad:

  individual line broadening parameter for broad (Lip / MM) signals. If
  NULL then defaults to lb_reg.

- broad_glb:

  apply global linewidth parameter to broad signals with names starting
  with Lip or MM.

## Value

full list of options.

## Examples

``` r
opts <- abfit_opts(ppm_left = 4.2, noise_region = c(-1, -3))
```
