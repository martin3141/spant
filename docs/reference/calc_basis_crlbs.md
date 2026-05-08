# Estimate the CRLB for each element in a basis set.

Estimate the CRLB for each element in a basis set.

## Usage

``` r
calc_basis_crlbs(
  basis,
  xlim = c(4, 0.2),
  zf = TRUE,
  sd = 1,
  bl_comp_pppm = NULL
)
```

## Arguments

- basis:

  basis_set object.

- xlim:

  spectral range to use in ppm.

- zf:

  zero-fill the basis set.

- sd:

  standard deviation of the noise.

- bl_comp_pppm:

  number spline baseline components to append per-ppm.

## Value

a vector of predicted errors.
