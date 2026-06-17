# Standard SVS 1H brain analysis pipeline.

Note this function is under active development and liable to changes.

## Usage

``` r
fit_svs(
  input,
  w_ref = NULL,
  output_dir = NULL,
  mri = NULL,
  mri_seg = NULL,
  deface = FALSE,
  segment_t1 = FALSE,
  segment_t1_method = "rpyants",
  external_basis = NULL,
  append_external_basis = FALSE,
  p_vols = NULL,
  format = NULL,
  pul_seq = NULL,
  TE = NULL,
  TR = NULL,
  TE1 = NULL,
  TE2 = NULL,
  TE3 = NULL,
  TM = NULL,
  append_basis = NULL,
  remove_basis = NULL,
  remove_external_basis = NULL,
  pre_align = TRUE,
  pre_align_max_shift = 40,
  pre_align_ref_freq = c(2.01, 3.03, 3.22),
  pre_align_ref_amp = 1,
  dfp_corr = TRUE,
  dfp_corr_ref_subset = NULL,
  output_ratio = NULL,
  ecc = NULL,
  hsvd_width = NULL,
  decimate = FALSE,
  trunc_fid_pts = NULL,
  fit_method = NULL,
  fit_opts = NULL,
  fit_subset = NULL,
  w_ref_subset = NULL,
  legacy_ws = FALSE,
  w_att = 0.7,
  w_conc = 35880,
  use_basis_cache = "auto",
  summary_measures = NULL,
  dyn_av_block_size = NULL,
  dyn_av_scheme = NULL,
  dyn_av_scheme_file = NULL,
  dyn_basis_lb = NULL,
  dyn_basis_lg = NULL,
  lcm_bin_path = NULL,
  plot_ppm_xlim = NULL,
  extra_output = FALSE,
  verbose = FALSE,
  return_fit = FALSE,
  write_preproc_metab_path = NULL,
  overwrite = FALSE
)
```

## Arguments

- input:

  path or mrs_data object containing MRS data.

- w_ref:

  path or mrs_data object containing MRS water reference data.

- output_dir:

  directory path to output fitting results.

- mri:

  filepath or nifti object containing anatomical MRI data.

- mri_seg:

  filepath or nifti object containing segmented MRI data.

- deface:

  option to apply faceoff to the mri input. Defaults to FALSE.

- segment_t1:

  segment the t1 weighted mri file with ANTs and use the results to
  perform partial volume correction. Defaults to FALSE.

- segment_t1_method:

  one of : "rpyants" (default), "ants" or "fslr".

- external_basis:

  precompiled basis set object to use for analysis.

- append_external_basis:

  append the external basis with the internally generated one. Useful
  for adding experimentally acquired baseline signals to internally
  simulated basis sets. Defaults to FALSE - meaning only signals from
  the external basis will be used in analysis.

- p_vols:

  a numeric vector of partial volumes expressed as percentages. Defaults
  to 100% white matter. A voxel containing 100% gray matter tissue would
  use : p_vols = c(WM = 0, GM = 100, CSF = 0).

- format:

  Override automatic data format detection. See format argument in
  [`read_mrs()`](https://martin3141.github.io/spant/reference/read_mrs.md)
  for permitted values.

- pul_seq:

  Pulse sequence to use for basis simulation. Can be one of the
  following values : "press", "press_ideal", "press_shaped", "steam" or
  "slaser". If "press" then "press_ideal" will be assumed unless the
  magnetic field is stronger that 2.8 Tesla, "press_shaped" will be
  assumed for 2.9 Tesla and above.

- TE:

  metabolite mrs data echo time in seconds. If not supplied this will be
  guessed from the metab data file.

- TR:

  metabolite mrs data repetition time in seconds. If not supplied this
  will be guessed from the metab data file.

- TE1:

  PRESS or sLASER sequence timing parameter in seconds.

- TE2:

  PRESS or sLASER sequence timing parameter in seconds.

- TE3:

  sLASER sequence timing parameter in seconds.

- TM:

  STEAM mixing time parameter in seconds.

- append_basis:

  names of extra signals to add to the default basis. Eg append_basis =
  c("peth", "cit"). Use get_mol_names() function to print all available
  signals. Cannot be used with precompiled basis sets.

- remove_basis:

  grep expression to match names of signals to remove from the basis.
  For example: use "^lac\$\|^ala\$" to remove lactate and alanine; "\*"
  to remove all signals and "^mm\|^lip" to remove all macromolecular and
  lipid signals. This operation is performed before signals are added
  with append_basis. Cannot be used with precompiled/exernal basis sets.

- remove_external_basis:

  grep expression to match names of signals to remove from the external
  basis. For example: use "^Lac\$\|^Ala\$" to remove lactateand alanine
  and "^MM\|^Lip" to remove all macromolecular and lipid signals.

- pre_align:

  perform simple frequency alignment to known reference peaks.

- pre_align_max_shift:

  maximum allowable shift in Hz. Defaults to 40 Hz.

- pre_align_ref_freq:

  reference frequency in ppm units. More than one frequency may be
  specified. Defaults to : c(2.01, 3.03, 3.22).

- pre_align_ref_amp:

  reference frequency relative amplitudes. Defaults to 1, corresponding
  to equal amplitudes.

- dfp_corr:

  perform dynamic frequency and phase correction using the RATS method.

- dfp_corr_ref_subset:

  specify a subset of dynamics to use as reference scans for dynamic
  frequency and phase correction. For example, 1:16 would use the mean
  of the first 16 dynamic scans as a reference spectrum.

- output_ratio:

  optional string to specify a metabolite ratio to output. Defaults to
  "tCr". Multiple metabolites may be specified for multiple outputs. Set
  to NA to omit.

- ecc:

  option to perform water reference based eddy current correction,
  default is to not apply unless the is GE format.

- hsvd_width:

  set the width of the HSVD filter in Hz. Note the applied width is
  between -width and +width Hz, with 0 Hz being defined at the centre of
  the spectral width. Default is disabled (set to NULL), 30 Hz is a
  reasonable value.

- decimate:

  option on decimate the data by a factor of 2 before analysis. Defaults
  to FALSE.

- trunc_fid_pts:

  number of points to truncate the input data by in the time-domain.
  E.g. setting to 1024 will ensure data with more time-domain points
  will be truncated to a length of 1024. Defaults to NULL, where
  truncation is not performed.

- fit_method:

  can be "ABFIT-REG" or "LCMODEL. Defaults to "ABFIT-REG".

- fit_opts:

  options to pass to the fitting method.

- fit_subset:

  specify a subset of dynamics to analyse, for example 1:16 would only
  fit the first 16 dynamic scans.

- w_ref_subset:

  specify a subset of dynamics to use as water reference following rats
  alignment and averaging. For example c(1, 2) would only use the
  average of the first two dynamic water reference scans. Default value
  of NULL uses all available water reference scans.

- legacy_ws:

  perform and output legacy water scaling compatible with default
  LCModel and TARQUIN behaviour. See w_att and w_conc arguments to
  change the default assumptions. Default value is FALSE.

- w_att:

  water attenuation factor (default = 0.7) for legacy water scaling.
  Assumes water T2 of 80ms and a TE = 30 ms. exp(-30ms / 80ms) ~ 0.7.

- w_conc:

  assumed water concentration (default = 35880) for legacy water
  scaling. Default value corresponds to typical white matter. Set to
  43300 for gray matter, and 55556 for phantom measurements.

- use_basis_cache:

  Pre-cache basis sets to reduce analysis speed. Can be one of the
  following : "auto", "all" or "none". The default value of "auto" will
  only use the cache for 3T PRESS - which generally requires more
  detailed simulation due to high CSD.

- summary_measures:

  output an additional table with a subset of metabolite levels, eg
  c("tNAA", "tNAA/tCr", "tNAA/tCho", "Lac/tNAA").

- dyn_av_block_size:

  perform temporal averaging with the specified block size. Defaults to
  NULL, eg average across all dynamic scans.

- dyn_av_scheme:

  a numeric vector of sequential integers (starting at 1), with the same
  length as the number of dynamic scans in the metabolite data. For
  example: c(1, 1, 2, 1, 1, 3, 1, 1).

- dyn_av_scheme_file:

  a file path containing a single column of sequential integers
  (starting at 1) with the same length as the number of dynamic scans in
  the metabolite data. File may be formatted as .xlsx, .xls, text or csv
  format.

- dyn_basis_lb:

  dynamic basis line-broadening to apply in Hz.

- dyn_basis_lg:

  dynamic basis Lorentz-Gauss lineshape factor between 0 and 1. Defaults
  to 0, pure Lorentzian.

- lcm_bin_path:

  set the path to LCModel binary.

- plot_ppm_xlim:

  plotting ppm axis limits in the html results. results.

- extra_output:

  write extra output files for generating custom plots. Defaults to
  FALSE.

- verbose:

  output potentially useful information.

- return_fit:

  return a fit object, defaults to FALSE.

- write_preproc_metab_path:

  path to write the preprocessed metabolite data in NIfTI format.

- overwrite:

  overwrite existing fitting result files, defaults to FALSE.

## Examples

``` r
metab <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
                     package = "spant")
w_ref <- system.file("extdata", "philips_spar_sdat_W.SDAT",
                     package = "spant")
out_dir <- file.path("~", "fit_svs_result")
if (FALSE) { # \dontrun{
fit_result <- fit_svs(metab, w_ref, out_dir)
} # }
```
