# Standard SVS 1H brain analysis pipeline.

Note this function is still under development and liable to changes.

## Usage

``` r
svs_1h_brain_analysis_dev(
  metab,
  w_ref = NULL,
  output_dir = NULL,
  basis = NULL,
  p_vols = NULL,
  append_basis = NULL,
  remove_basis = NULL,
  dfp_corr = FALSE,
  omit_bad_dynamics = FALSE,
  te = NULL,
  tr = NULL,
  output_ratio = "tCr",
  ecc = FALSE,
  abfit_opts = NULL,
  verbose = FALSE
)
```

## Arguments

- metab:

  filepath or mrs_data object containing MRS metabolite data.

- w_ref:

  filepath or mrs_data object containing MRS water reference data.

- output_dir:

  directory path to output fitting results.

- basis:

  precompiled basis set object to use for analysis.

- p_vols:

  a numeric vector of partial volumes expressed as percentages. Defaults
  to 100% white matter. A voxel containing 100% gray matter tissue would
  use : p_vols = c(WM = 0, GM = 100, CSF = 0).

- append_basis:

  names of extra signals to add to the default basis. Eg append_basis =
  c("peth", "cit"). Cannot be used with precompiled basis sets.

- remove_basis:

  names of signals to remove from the basis. Cannot be used with
  precompiled basis sets.

- dfp_corr:

  perform dynamic frequency and phase correction using the RATS method.

- omit_bad_dynamics:

  detect and remove bad dynamics.

- te:

  metabolite mrs data echo time in seconds. If not supplied this will be
  guessed from the metab data file.

- tr:

  metabolite mrs data repetition time in seconds. If not supplied this
  will be guessed from the metab data file.

- output_ratio:

  optional string to specify a metabolite ratio to output. Defaults to
  "tCr" and multiple metabolites may be specified for multiple outputs.
  Set as NULL to omit.

- ecc:

  option to perform water reference based eddy current correction,
  defaults to FALSE.

- abfit_opts:

  options to pass to ABfit.

- verbose:

  output potentially useful information.

## Examples

``` r
metab <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
                     package = "spant")
w_ref <- system.file("extdata", "philips_spar_sdat_W.SDAT",
                     package = "spant")
if (FALSE) { # \dontrun{
fit_result <- svs_1h_brain_analysis(metab, w_ref, "fit_res_dir")
} # }
```
