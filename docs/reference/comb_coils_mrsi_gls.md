# Combine MRSI coil data using the GLS method presented by An et al JMRI 37:1445-1450 (2013).

Combine MRSI coil data using the GLS method presented by An et al JMRI
37:1445-1450 (2013).

## Usage

``` r
comb_coils_mrsi_gls(
  metab,
  noise_pts = 256,
  noise_mrs = NULL,
  bc_poly_noise = 2
)
```

## Arguments

- metab:

  MRSI data containing metabolite data.

- noise_pts:

  number of points from the end of the FIDs to use for noise covariance
  estimation.

- noise_mrs:

  MRS data containing noise information for each coil.

- bc_poly_noise:

  baseline correct the noise samples with a polynomial in the
  time-domain. Defaults to 2, which performs a second-order polynomial
  correction. Set to NULL to disable.

## Value

coil combined MRSI data.
