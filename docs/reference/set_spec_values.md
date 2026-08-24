# Set spectral data points to a given value.

Set spectral data points to a given value.

## Usage

``` r
set_spec_values(mrs_data, xlim = NULL, value = 0, invert = FALSE)
```

## Arguments

- mrs_data:

  data to be processed.

- xlim:

  spectral region (in ppm) of data points to set. If NULL (default) all
  spectral points will be set.

- value:

  value to set the spectral region. Defaults to zero.

- invert:

  set all spectral points *outside* the spectral region to the given
  value.

## Value

mrs_data object with modified spectral data points.
