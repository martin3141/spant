# Perform zeroth-order phase correction based on the minimisation of the squared difference between the real and magnitude components of the spectrum.

Perform zeroth-order phase correction based on the minimisation of the
squared difference between the real and magnitude components of the
spectrum.

## Usage

``` r
auto_phase(mrs_data, xlim = c(4, 1.8), smo_ppm_sd = 0.05, ret_phase = FALSE)
```

## Arguments

- mrs_data:

  an object of class `mrs_data`.

- xlim:

  frequency range (default units of PPM) to including in the phase.

- smo_ppm_sd:

  Gaussian smoother sd in ppm units.

- ret_phase:

  return phase values (logical).

## Value

MRS data object and phase values (optional).
