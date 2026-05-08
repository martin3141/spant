# Fit a single asymmetric pseudo-Voigt resonance in the frequency domain.

Fit a single asymmetric pseudo-Voigt resonance in the frequency domain.

## Usage

``` r
fit_asy_pvoigt(
  mrs_data,
  freq_ppm = 4.65,
  xlim = c(5.2, 4.1),
  lg_limits = c(0, 1)
)
```

## Arguments

- mrs_data:

  data containing the resonance to be fit.

- freq_ppm:

  frequency estimate (in ppm) for the resonance to be fitted.

- xlim:

  spectral range (in ppm) where the fit will be evaluated.

- lg_limits:

  lg lineshape parameter limits.

## Value

list of fitting results and parameters.
