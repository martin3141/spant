# Corrected zero order phase and chemical shift offset in 1H MRS data from the brain.

Corrected zero order phase and chemical shift offset in 1H MRS data from
the brain.

## Usage

``` r
phase_ref_1h_brain(
  mrs_data,
  mean_ref = FALSE,
  ret_corr_only = TRUE,
  xlim = c(4, 1.9),
  p_deg = 3,
  sp_N = 2,
  basis_type = "poly"
)
```

## Arguments

- mrs_data:

  MRS data to be corrected.

- mean_ref:

  apply the phase and offset of the mean spectrum to all others. Default
  is FALSE.

- ret_corr_only:

  return the corrected data only.

- xlim:

  frequency range in ppm to consider.

- p_deg:

  polynomial baseline order.

- sp_N:

  number of spline functions, note the true number will be sp_N +
  sp_deg.

- basis_type:

  may be one of "poly" or "spline".

## Value

corrected MRS data.
