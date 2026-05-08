# Preprocess and perform quality assessment of one or more SVS data sets.

Preprocess and perform quality assessment of one or more SVS data sets.

## Usage

``` r
preproc_svs_dataset(
  paths,
  labels = NULL,
  output_dir = "spant_analysis",
  exclude_labels = NULL,
  overwrite = FALSE,
  ref_inds = NULL,
  dyn_water_ref_paths = NULL,
  Ntrans = NULL,
  TR = NULL,
  return_results = FALSE
)
```

## Arguments

- paths:

  paths to the MRS data file or IMA directory.

- labels:

  labels to describe each data set.

- output_dir:

  output directory.

- exclude_labels:

  vector of labels of scans to exclude, eg poor quality data.

- overwrite:

  overwrite saved results, defaults to FALSE.

- ref_inds:

  a vector of 1-based indices for any water reference dynamic scans.

- dyn_water_ref_paths:

  paths to interleaved water reference scans.

- Ntrans:

  override the number of transients set in the mrs_data file.

- TR:

  override the TR (in seconds) set in the mrs_data file.

- return_results:

  function will return key outputs, defaults to FALSE.
