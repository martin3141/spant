# Combine all fitting data points into a single data frame.

Combine all fitting data points into a single data frame.

## Usage

``` r
comb_fit_tables(fit_res, inc_basis_sigs = FALSE, inc_indices = TRUE)
```

## Arguments

- fit_res:

  a single fit_result object.

- inc_basis_sigs:

  include the individual fitting basis signals in the output table,
  defaults to FALSE.

- inc_indices:

  include indices such as X, Y and coil in the output, defaults to TRUE.
  These are generally not useful for SVS analysis.

## Value

a data frame containing the fit data points.
