# Combine MRSI coil data using the GLS method presented by An et al JMRI 37:1445-1450 (2013).

Combine MRSI coil data using the GLS method presented by An et al JMRI
37:1445-1450 (2013).

## Usage

``` r
comb_coils_mrsi_gls(metab, noise_pts = 30, noise_mrs = NULL)
```

## Arguments

- metab:

  MRSI data containing metabolite data.

- noise_pts:

  number of points from the end of the FIDs to use for noise covariance
  estimation.

- noise_mrs:

  MRS data containing noise information for each coil.

## Value

coil combined MRSI data.
