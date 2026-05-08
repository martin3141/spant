# Collapse MRS data by concatenating spectra along the dynamic dimension.

Collapse MRS data by concatenating spectra along the dynamic dimension.

## Usage

``` r
collapse_to_dyns(x, rm_masked = FALSE)

# S3 method for class 'mrs_data'
collapse_to_dyns(x, rm_masked = FALSE)

# S3 method for class 'fit_result'
collapse_to_dyns(x, rm_masked = FALSE)
```

## Arguments

- x:

  data object to be collapsed (mrs_data or fit_result object).

- rm_masked:

  remove masked dynamics from the output.

## Value

collapsed data with spectra or fits concatenated along the dynamic
dimension.
