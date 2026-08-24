# Set all spectral data points below a given percentage of the maximum to zero.

Set all spectral data points below a given percentage of the maximum to
zero.

## Usage

``` r
zero_spec_threshold(mrs_data, percent_max)
```

## Arguments

- mrs_data:

  data to be processed.

- percent_max:

  spectral data points smaller than this percentage value will be set to
  zero.

## Value

mrs_data object with zeroed spectral data points.
