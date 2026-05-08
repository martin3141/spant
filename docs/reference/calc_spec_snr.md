# Calculate the spectral SNR.

SNR is defined as the maximum signal value divided by the standard
deviation of the noise.

## Usage

``` r
calc_spec_snr(
  mrs_data,
  sig_region = c(4, 0.5),
  noise_region = c(-0.5, -2.5),
  p_order = 2,
  interp_f = 4,
  full_output = FALSE
)
```

## Arguments

- mrs_data:

  an object of class `mrs_data`.

- sig_region:

  a ppm region to define where the maximum signal value should be
  estimated.

- noise_region:

  a ppm region to defined where the noise level should be estimated.

- p_order:

  polynomial order to fit to the noise region before estimating the
  standard deviation.

- interp_f:

  interpolation factor to improve detection of the highest signal value.

- full_output:

  output signal, noise and SNR values separately.

## Value

an array of SNR values.

## Details

The mean noise value is subtracted from the maximum signal value to
reduce DC offset bias. A polynomial detrending fit (second order by
default) is applied to the noise region before the noise standard
deviation is estimated.
