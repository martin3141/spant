# Apply Voigt line-broadening to match a reference spectrum.

Apply Voigt line-broadening to match a reference spectrum.

## Usage

``` r
match_lineshape(
  mrs_data,
  ref,
  xlim,
  init_lb = 0.2,
  init_lg = 0.5,
  init_amp = 1,
  min_lb = 0,
  max_lb = Inf,
  min_lg = 0,
  max_lg = 1,
  min_amp = 0.1,
  max_amp = 2,
  amp_optim = TRUE
)
```

## Arguments

- mrs_data:

  data to be broadened, note the linewidth of this spectrum must be
  narrower than the ref spectrum.

- ref:

  reference data to match.

- xlim:

  spectral region to match, eg c(5.2, 4.1) could be used to match two
  water resonances.

- init_lb:

  initial value for the amount of line-broadening to apply (Hz).

- init_lg:

  initial value for the Lorentz-Gauss lineshape parameter.

- init_amp:

  initial value for the amplitude parameter.

- min_lb:

  minimum value for the amount of line-broadening to apply (Hz).

- max_lb:

  maximum value for the amount of line-broadening to apply (Hz).

- min_lg:

  minimum value for the Lorentz-Gauss lineshape parameter.

- max_lg:

  maximum value for the Lorentz-Gauss lineshape parameter.

- min_amp:

  minimum value for the amplitude parameter.

- max_amp:

  maximum value for the amplitude parameter.

- amp_optim:

  flag to include amplitude adjustment in the optimisation procedure.
  Defaults to TRUE.

## Value

a list containing the matched mrs_data, difference spectra and
optimisation results.
