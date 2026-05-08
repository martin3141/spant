# Perform l2 regularisation artefact suppression.

Perform l2 regularisation artefact suppression using the method proposed
by Bilgic et al. JMRI 40(1):181-91 2014.

## Usage

``` r
l2_reg(
  mrs_data,
  thresh = 0.05,
  b = 1e-11,
  A = NA,
  xlim = NA,
  thresh_xlim = NULL,
  A_append = NULL,
  ret_norms = FALSE
)
```

## Arguments

- mrs_data:

  input data for artefact suppression.

- thresh:

  threshold parameter to extract lipid signals from mrs_data based on
  the spectral integration of the thresh_xlim region in magnitude mode.

- b:

  regularisation parameter.

- A:

  set of spectra containing the artefact basis signals. The thresh
  parameter is ignored when A is specified.

- xlim:

  spectral limits in ppm to restrict the reconstruction range. Defaults
  to the full spectral width.

- thresh_xlim:

  spectral limits in ppm to integrate for the threshold map.

- A_append:

  additional spectra to append to the A basis.

- ret_norms:

  return the residual norm and solution norms.

## Value

l2 reconstructed mrs_data object.
