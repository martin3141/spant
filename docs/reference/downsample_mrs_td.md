# Downsample an MRS signal by a factor of 2 by removing every other data point in the time-domain. Note, signals outside the new sampling frequency will be aliased.

Downsample an MRS signal by a factor of 2 by removing every other data
point in the time-domain. Note, signals outside the new sampling
frequency will be aliased.

## Usage

``` r
downsample_mrs_td(mrs_data)
```

## Arguments

- mrs_data:

  MRS data object.

## Value

downsampled data.
