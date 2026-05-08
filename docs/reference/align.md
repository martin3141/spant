# Align spectra to a reference frequency using a convolution based method.

Align spectra to a reference frequency using a convolution based method.

## Usage

``` r
align(
  mrs_data,
  ref_freq = 4.65,
  ref_amp = 1,
  zf_factor = 2,
  lb = 2,
  max_shift = 20,
  ret_df = FALSE,
  mean_dyns = FALSE
)
```

## Arguments

- mrs_data:

  data to be aligned.

- ref_freq:

  reference frequency in ppm units. More than one frequency may be
  specified.

- ref_amp:

  amplitude value for the reference signal. More than one value may be
  specified to match the number of ref_freq signals.

- zf_factor:

  zero filling factor to increase alignment resolution.

- lb:

  line broadening to apply to the reference signal.

- max_shift:

  maximum allowable shift in Hz.

- ret_df:

  return frequency shifts in addition to aligned data (logical).

- mean_dyns:

  align the mean spectrum and apply the same shift to each dynamic.

## Value

aligned data object.
