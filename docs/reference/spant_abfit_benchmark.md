# Simulate and fit some spectra with ABfit for benchmarking purposes. Basic timing and performance metrics will be printed.

Simulate and fit some spectra with ABfit for benchmarking purposes.
Basic timing and performance metrics will be printed.

## Usage

``` r
spant_abfit_benchmark(noise_reps = 10, return_res = FALSE, opts = abfit_opts())
```

## Arguments

- noise_reps:

  number of spectra to fit with differing noise samples.

- return_res:

  return a list of fit_result objects.

- opts:

  ABfit options structure.
