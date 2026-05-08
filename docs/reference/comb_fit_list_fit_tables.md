# Combine all fitting data points from a list of fits into a single data frame.

Combine all fitting data points from a list of fits into a single data
frame.

## Usage

``` r
comb_fit_list_fit_tables(
  fit_list,
  add_extra = TRUE,
  harmonise_ppm = TRUE,
  inc_basis_sigs = FALSE,
  inc_indices = TRUE,
  add_res_id = TRUE
)
```

## Arguments

- fit_list:

  list of fit_result objects.

- add_extra:

  add variables in the extra data frame to the output (TRUE).

- harmonise_ppm:

  ensure the ppm scale for each fit is identical to the first.

- inc_basis_sigs:

  include the individual fitting basis signals in the output table,
  defaults to FALSE.

- inc_indices:

  include indices such as X, Y and coil in the output, defaults to TRUE.
  These are generally not useful for SVS analysis.

- add_res_id:

  add a res_id column to the output to distinguish between datasets.

## Value

a data frame containing the fit data points.
