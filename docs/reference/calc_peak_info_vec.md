# Calculate the FWHM of a peak from a vector of intensity values.

Calculate the FWHM of a peak from a vector of intensity values.

## Usage

``` r
calc_peak_info_vec(data_pts, interp_f)
```

## Arguments

- data_pts:

  input vector.

- interp_f:

  interpolation factor to improve the FWHM estimate.

## Value

a vector of: x position of the highest data point, maximum peak value in
the y axis, FWHM in the units of data points.
