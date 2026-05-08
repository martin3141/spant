# Perform zeroth-order phase correction based on expected baseline regions.

Default arguments are appropriate for water suppressed 1H MRS of the
brain. The following options are suitable for a water signal at 4.65
ppm: ppm_start = c(5.2, 4.0), ppm_end = c(5.3, 4.1), xlim = c(5.2, 4.1)

## Usage

``` r
auto_phase_bl(
  mrs_data,
  ppm_start = c(0.5, 4),
  ppm_end = c(0, 4.2),
  xlim = c(1.8, 4),
  mean_dyns = FALSE,
  ret_phase = FALSE
)
```

## Arguments

- mrs_data:

  an object of class `mrs_data`.

- ppm_start:

  a vector of ppm values designating baseline regions.

- ppm_end:

  a vectors of ppm values designating baseline regions.

- xlim:

  region containing signal of interest, eg strong metabolite resonances.

- mean_dyns:

  phase the mean spectrum and apply the same value to each dynamic.

- ret_phase:

  return phase values (logical).

## Value

MRS data object with corrected phase.
