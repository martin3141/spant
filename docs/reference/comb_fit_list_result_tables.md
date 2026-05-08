# Combine the fit result tables from a list of fit results.

Combine the fit result tables from a list of fit results.

## Usage

``` r
comb_fit_list_result_tables(fit_list, add_extra = TRUE, add_res_id = TRUE)
```

## Arguments

- fit_list:

  a list of fit_result objects.

- add_extra:

  add variables in the extra data frame to the output (TRUE).

- add_res_id:

  add a res_id column to the output to distinguish between datasets.

## Value

a data frame combine all fit result tables with an additional id column
to differentiate between data sets. Any variables in the extra data
frame may be optionally added to the result.
