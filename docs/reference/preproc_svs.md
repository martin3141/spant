# Preprocess and perform quality assessment of a single SVS data set.

Preprocess and perform quality assessment of a single SVS data set.

## Usage

``` r
preproc_svs(
  path,
  label = NULL,
  output_dir = NULL,
  ref_inds = NULL,
  dyn_water_ref_path = NULL,
  Ntrans = NULL,
  TR = NULL
)
```

## Arguments

- path:

  path to the MRS data file or IMA directory.

- label:

  a label to describe the data set.

- output_dir:

  output directory.

- ref_inds:

  a vector of 1-based indices for any water reference dynamic scans.

- dyn_water_ref_path:

  path to interleaved water reference data file or IMA directory.

- Ntrans:

  override the number of transients set in the mrs_data file.

- TR:

  override the TR (in seconds) set in the mrs_data file.
