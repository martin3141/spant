# Fit a T2 relaxation curve, from multiple TEs, to a set of amplitudes.

Fit a T2 relaxation curve, from multiple TEs, to a set of amplitudes.

## Usage

``` r
fit_t2_te_array(
  te_vec,
  amp_vec,
  lower = 0,
  upper = 10,
  output_fit_res = 0.01,
  ret_full = TRUE
)
```

## Arguments

- te_vec:

  vector of TE values in seconds.

- amp_vec:

  vector of amplitudes.

- lower:

  minimum allowable T2 value.

- upper:

  maximum allowable T2 value.

- output_fit_res:

  temporal resolution (seconds) of the ideal output relaxation curve.

- ret_full:

  return full fitting information including ideal relaxation curve.

## Value

a list containing relaxation parameters and an ideal curve for fit
evaluation.
