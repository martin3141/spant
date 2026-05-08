# Scale an mrs_data object by a scalar or vector or amplitudes.

Scale an mrs_data object by a scalar or vector or amplitudes.

## Usage

``` r
scale_mrs_amp(mrs_data, amp)
```

## Arguments

- mrs_data:

  data to be scaled.

- amp:

  multiplicative factor, must have length equal to 1 or Nspec(mrs_data).

## Value

mrs_data object multiplied by the amplitude scale factor.
