# Fit and subtract a smoothing spline to each spectrum in a dataset.

Fit and subtract a smoothing spline to each spectrum in a dataset.

## Usage

``` r
bc_spline(mrs_data, spar = 0.5, nknots = 100)
```

## Arguments

- mrs_data:

  mrs_data object.

- spar:

  smoothing parameter typically between 0 and 1.

- nknots:

  number of spline knots.

## Value

smoothing spline subtracted data.
