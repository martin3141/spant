# Estimate the standard deviation of the noise from a segment of an mrs_data object.

Estimate the standard deviation of the noise from a segment of an
mrs_data object.

## Usage

``` r
est_noise_sd(mrs_data, n = 100, offset = 100, p_order = 2)
```

## Arguments

- mrs_data:

  MRS data object.

- n:

  number of data points (taken from the end of array) to use in the
  estimation.

- offset:

  number of final points to exclude from the calculation.

- p_order:

  polynomial order to fit to the data before estimating the standard
  deviation.

## Value

standard deviation array.
