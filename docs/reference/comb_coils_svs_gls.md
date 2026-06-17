# Combine SVS coil data using the GLS method presented by An et al JMRI 37:1445-1450 (2013).

Combine SVS coil data using the GLS method presented by An et al JMRI
37:1445-1450 (2013).

## Usage

``` r
comb_coils_svs_gls(
  metab,
  ref = NULL,
  noise_pts = 256,
  noise_mrs = NULL,
  use_mean_sens = TRUE,
  use_ref_sens = FALSE,
  change_order = NULL,
  no_noise_corr = FALSE,
  bc_poly_noise = 2
)
```

## Arguments

- metab:

  MRS data containing metabolite data.

- ref:

  MRS data containing reference data (optional).

- noise_pts:

  number of points from the end of the FIDs to use for noise covariance
  estimation.

- noise_mrs:

  MRS data containing noise information for each coil.

- use_mean_sens:

  use the dynamic mean to estimate coil sensitivities.

- use_ref_sens:

  use the reference data to estimate coil sensitivities.

- change_order:

  change the order of the noise matrix reshaping, for testing purposes
  only.

- no_noise_corr:

  assume noise is uncorrelated between coils, for testing purposes only.
  Defaults to FALSE.

- bc_poly_noise:

  baseline correct the noise samples with a polynomial in the
  time-domain. Defaults to 2, which performs a second-order polynomial
  correction. Set to NULL to disable.

## Value

coil combined MRS data.
