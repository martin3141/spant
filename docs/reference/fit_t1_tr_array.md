# Fit a T1 recovery curve, from multiple TRs, to a set of amplitudes.

Fit a T1 recovery curve, from multiple TRs, to a set of amplitudes.

## Usage

``` r
fit_t1_tr_array(
  tr_vec,
  amp_vec,
  lower = 0,
  upper = 10,
  output_fit_res = 0.01,
  ret_full = TRUE
)
```

## Arguments

- tr_vec:

  vector of TR values in seconds.

- amp_vec:

  vector of amplitudes.

- lower:

  minimum allowable T1 value.

- upper:

  maximum allowable T1 value.

- output_fit_res:

  temporal resolution (seconds) of the ideal output relaxation curve.

- ret_full:

  return full fitting information including ideal relaxation curve.

## Value

a list containing relaxation parameters and an ideal curve for fit
evaluation.
