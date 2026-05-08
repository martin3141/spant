# Standard SVS 1H brain analysis pipeline.

Standard SVS 1H brain analysis pipeline.

## Usage

``` r
svs_1h_brain_analysis(
  metab,
  basis = NULL,
  w_ref = NULL,
  mri_seg = NULL,
  mri = NULL,
  output_dir = NULL,
  extra = NULL,
  decimate = NULL,
  rats_corr = TRUE,
  ecc = FALSE,
  comb_dyns = TRUE,
  hsvd_filt = FALSE,
  scale_amps = TRUE,
  te = NULL,
  tr = NULL,
  preproc_only = FALSE,
  method = "ABFIT",
  opts = NULL
)
```

## Arguments

- metab:

  filepath or mrs_data object containing MRS metabolite data.

- basis:

  basis set object to use for analysis.

- w_ref:

  filepath or mrs_data object containing MRS water reference data.

- mri_seg:

  filepath or nifti object containing segmented MRI data.

- mri:

  filepath or nifti object containing anatomical MRI data.

- output_dir:

  directory path to output fitting results.

- extra:

  data.frame with one row containing additional information to be
  attached to the fit results table.

- decimate:

  option to decimate the input data by a factor of two. The default
  value of NULL does not perform decimation unless the spectral width is
  greater than 20 PPM.

- rats_corr:

  option to perform rats correction, defaults to TRUE.

- ecc:

  option to perform water reference based eddy current correction,
  defaults to FALSE.

- comb_dyns:

  option to combine dynamic scans, defaults to TRUE.

- hsvd_filt:

  option to apply hsvd water removal, defaults to FALSE.

- scale_amps:

  option to scale metabolite amplitude estimates, defaults to TRUE.

- te:

  metabolite mrs data echo time in seconds.

- tr:

  metabolite mrs data repetition time in seconds.

- preproc_only:

  only perform the preprocessing steps and omit fitting. The
  preprocessed metabolite data will be returned in this case.

- method:

  analysis method to use, see fit_mrs help.

- opts:

  options to pass to the analysis method.

## Value

a fit_result or mrs_data object depending on the preproc_only option.
