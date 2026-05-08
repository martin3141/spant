# Perform a zeroth order phase correction based on the phase of the first data point in the time-domain.

Perform a zeroth order phase correction based on the phase of the first
data point in the time-domain.

## Usage

``` r
fp_phase_correct(mrs_data, ret_phase = FALSE)
```

## Arguments

- mrs_data:

  MRS data to be corrected.

- ret_phase:

  return phase values (logical).

## Value

corrected data or a list with corrected data and optional phase values.
