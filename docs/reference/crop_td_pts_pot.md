# Crop `mrs_data` object data points in the time-domain rounding down to the next smallest power of two (pot). Data that already has a pot length will not be changed.

Crop `mrs_data` object data points in the time-domain rounding down to
the next smallest power of two (pot). Data that already has a pot length
will not be changed.

## Usage

``` r
crop_td_pts_pot(mrs_data)
```

## Arguments

- mrs_data:

  MRS data.

## Value

cropped `mrs_data` object.
