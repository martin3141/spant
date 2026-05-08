# Edited SVS 1H brain analysis pipeline.

Note this function is still under development and liable to changes.

## Usage

``` r
fit_svs_edited(
  input,
  w_ref = NULL,
  output_dir = NULL,
  mri = NULL,
  mri_seg = NULL,
  deface = FALSE,
  segment_t1 = FALSE,
  segment_t1_method = "ants",
  external_basis = NULL,
  p_vols = NULL,
  format = NULL,
  editing_type = "gaba_1.9",
  editing_scheme = NULL,
  invert_edit_on = NULL,
  invert_edit_off = NULL,
  pul_seq = NULL,
  TE = NULL,
  TR = NULL,
  TE1 = NULL,
  TE2 = NULL,
  TE3 = NULL,
  TM = NULL,
  append_basis_ed_off = NULL,
  remove_basis_ed_off = NULL,
  pre_align = TRUE,
  dfp_corr = TRUE,
  output_ratio = NULL,
  ecc = FALSE,
  hsvd_width = NULL,
  decimate = FALSE,
  trunc_fid_pts = NULL,
  fit_opts_edited = NULL,
  fit_opts_ed_off = NULL,
  fit_subset = NULL,
  legacy_ws = FALSE,
  w_att = 0.7,
  w_conc = 35880,
  use_basis_cache = "auto",
  summary_measures = NULL,
  dyn_av_block_size = NULL,
  dyn_av_scheme = NULL,
  dyn_av_scheme_file = NULL,
  plot_ppm_xlim = NULL,
  extra_output = FALSE,
  verbose = FALSE,
  return_fit = FALSE,
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

  one of : "ants" (default), "rpyants" or "fslr".

- external_basis:

  precompiled basis set object to use for analysis.

- p_vols:

  a numeric vector of partial volumes expressed as percentages. Defaults
  to 100% white matter. A voxel containing 100% gray matter tissue would
  use : p_vols = c(WM = 0, GM = 100, CSF = 0).

- format:

  Override automatic data format detection. See format argument in
  [`read_mrs()`](https://martin3141.github.io/spant/reference/read_mrs.md)
  for permitted values.

- editing_type:

  can be one of : "gaba_1.9" or "gsh_4.54". Defaults to "gaba_1.9".

- editing_scheme:

  describes the dynamic data ordering. Can be one of: 'on-off-blocks',
  'on-off-interleaved', 'off-on-blocks' or 'off-on-interleaved'.

- invert_edit_on:

  set to TRUE to invert the edit-on sub-spectra.

- invert_edit_off:

  set to TRUE to invert the edit-off sub-spectra.

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

- append_basis_ed_off:

  names of extra signals to add to the default basis. Eg append_basis =
  c("peth", "cit"). Use get_mol_names() function to print all available
  signals. Cannot be used with precompiled basis sets.

- remove_basis_ed_off:

  grep expression to match names of signals to remove from the basis.
  For example: use "lac\|ala" to remove lactate and alanine; "\*" to
  remove all signals and "^mm\|^lip" to remove all macromolecular and
  lipid signals. This operation is performed before signals are added
  with append_basis_ed_off. Cannot be used with precompiled basis sets.

- pre_align:

  perform simple frequency alignment to known reference peaks.

- dfp_corr:

  perform dynamic frequency and phase correction using the RATS method.

- output_ratio:

  optional string to specify a metabolite ratio to output. Defaults to
  "tCr". Multiple metabolites may be specified for multiple outputs. Set
  to NA to omit.

- ecc:

  option to perform water reference based eddy current correction,
  defaults to FALSE.

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

- fit_opts_edited:

  options to pass to the fitting method for the edited spectrum.

- fit_opts_ed_off:

  options to pass to the fitting method for the edit-off spectrum.

- fit_subset:

  specify a subset of dynamics to analyse, for example 1:16 would only
  fit the first 16 dynamic scans.

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

- plot_ppm_xlim:

  plotting ppm axis limits in the html results. results.

- extra_output:

  write extra output files for generating custom plots. Defaults to
  FALSE.

- verbose:

  output potentially useful information.

- return_fit:

  return a fit object, defaults to FALSE.

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
