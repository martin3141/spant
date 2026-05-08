# Convert mrs_data object to a matrix, with spectral points in the column dimension and dynamics in the row dimension.

Convert mrs_data object to a matrix, with spectral points in the column
dimension and dynamics in the row dimension.

## Usage

``` r
mrs_data2mat(mrs_data, collapse = TRUE)
```

## Arguments

- mrs_data:

  MRS data object or list of MRS data objects.

- collapse:

  collapse all other dimensions along the dynamic dimension, eg a 16x16
  MRSI grid would be first collapsed across 256 dynamic scans.

## Value

MRS data matrix.
