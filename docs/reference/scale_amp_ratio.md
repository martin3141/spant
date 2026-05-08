# Scale fitted amplitudes to a ratio of signal amplitude.

Scale fitted amplitudes to a ratio of signal amplitude.

## Usage

``` r
scale_amp_ratio(fit_result, name, use_mean_value = FALSE)
```

## Arguments

- fit_result:

  a result object generated from fitting.

- name:

  the signal name to use as a denominator (usually, "tCr" or "tNAA").

- use_mean_value:

  scales the result by the mean of the signal when set to TRUE.

## Value

a `fit_result` object with a rescaled results table.
