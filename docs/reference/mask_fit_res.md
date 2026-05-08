# Mask fit result spectra depending on a vector of bool values.

Mask fit result spectra depending on a vector of bool values.

## Usage

``` r
mask_fit_res(fit_result, mask_vec, amps_only = FALSE)
```

## Arguments

- fit_result:

  fit result object to be masked.

- mask_vec:

  a Boolean vector with the same number of rows as there are rows in the
  results table.

- amps_only:

  only mask the amplitude and associated error estimate columns.

## Value

a masked fit result object.
