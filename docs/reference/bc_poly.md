# Fit and subtract a polynomial to each spectrum in a dataset.

Fit and subtract a polynomial to each spectrum in a dataset.

## Usage

``` r
bc_poly(mrs_data, p_deg = 1, fd = TRUE)
```

## Arguments

- mrs_data:

  mrs_data object.

- p_deg:

  polynomial degree.

- fd:

  perform correction in the frequency domain if TRUE, or time-domain if
  FALSE. Default is to perform in the frequency domain.

## Value

polynomial subtracted data.
