# Correct linear frequency drift.

Correct linear frequency drift.

## Usage

``` r
lofdc(
  mrs_data,
  max_hz_s = 0.1,
  tr = NULL,
  ret_corr_only = TRUE,
  outlier_thresh = 3,
  xlim = c(4, 0.5),
  order = 1
)
```

## Arguments

- mrs_data:

  MRS data to be corrected.

- max_hz_s:

  the maximum drift rate to search over.

- tr:

  mrs_data repetition time.

- ret_corr_only:

  return the corrected mrs_data object only.

- outlier_thresh:

  threshold to remove outliers.

- xlim:

  spectral width (in ppm) to evaluate outliers.

- order:

  correction order.

## Value

drift corrected mrs_data object.
