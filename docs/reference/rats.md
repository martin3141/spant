# Robust Alignment to a Target Spectrum (RATS).

Robust Alignment to a Target Spectrum (RATS).

## Usage

``` r
rats(
  mrs_data,
  ref = NULL,
  xlim = c(4, 0.5),
  max_shift = 20,
  p_deg = 2,
  sp_N = 2,
  sp_deg = 3,
  max_t = 0.2,
  basis_type = "poly",
  rescale_output = TRUE,
  phase_corr = TRUE,
  ret_corr_only = TRUE,
  zero_freq_shift_t0 = FALSE,
  list_mean_ref = TRUE,
  remove_freq_outliers = FALSE,
  freq_outlier_thresh = 3,
  remove_phase_outliers = FALSE,
  phase_outlier_thresh = 3,
  remove_amp_outliers = FALSE,
  amp_outlier_thresh = 3
)
```

## Arguments

- mrs_data:

  MRS data to be corrected.

- ref:

  optional MRS data to use as a reference, the mean of all dynamics is
  used if this argument is not supplied.

- xlim:

  optional frequency range to perform optimisation, set to NULL to use
  the full range.

- max_shift:

  maximum allowable frequency shift in Hz.

- p_deg:

  polynomial degree used for baseline modelling. Negative values disable
  baseline modelling.

- sp_N:

  number of spline functions, note the true number will be sp_N +
  sp_deg.

- sp_deg:

  degree of spline functions.

- max_t:

  truncate the FID when longer than max_t to reduce time taken, set to
  NULL to use the entire FID.

- basis_type:

  may be one of "poly" or "spline".

- rescale_output:

  rescale the bl_matched_spec and bl output to improve consistency
  between dynamic scans.

- phase_corr:

  apply phase correction (in addition to frequency). TRUE by default.

- ret_corr_only:

  return the corrected mrs_data object only.

- zero_freq_shift_t0:

  perform a linear fit to the frequency shifts and set the (linearly
  modeled) shift to be 0 Hz for the first dynamic scan.

- list_mean_ref:

  is ref is not specified and a list is provided, use the list mean scan
  as reference. Otherwise use the mean of each list element as it's own
  reference.

- remove_freq_outliers:

  remove dynamics based on their frequency shift.

- freq_outlier_thresh:

  threshold to remove frequency outliers.

- remove_phase_outliers:

  remove dynamics based on their phase shift.

- phase_outlier_thresh:

  threshold to remove phase outliers.

- remove_amp_outliers:

  remove dynamics based on their amplitude change.

- amp_outlier_thresh:

  threshold to remove amplitude outliers.

## Value

a list containing the corrected data; phase and shift values in units of
degrees and Hz respectively.
