# Baseline correction using the ALS method.

Eilers P. H. C. and Boelens H. F. M. (2005) Baseline correction with
asymmetric least squares smoothing. Leiden Univ. Medical Centre Report.

## Usage

``` r
bc_als(mrs_data, lambda = 10000, p = 0.001, ret_bc_only = TRUE)
```

## Arguments

- mrs_data:

  mrs_data object.

- lambda:

  controls the baseline flexibility.

- p:

  controls the penalty for negative data points.

- ret_bc_only:

  return the baseline corrected data only. When FALSE the baseline
  estimate and input data will be returned.

## Value

baseline corrected data.
