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
  use_mean_sens = TRUE
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

## Value

coil combined MRS data.
