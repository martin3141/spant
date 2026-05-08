# Add noise to an mrs_data object.

Add noise to an mrs_data object.

## Usage

``` r
add_noise(mrs_data, sd = 0.1, fd = TRUE)
```

## Arguments

- mrs_data:

  data to add noise to.

- sd:

  standard deviation of the noise.

- fd:

  generate the noise samples in the frequency-domain (TRUE) or
  time-domain (FALSE). This is required since the absolute value of the
  standard deviation of noise samples changes when data is Fourier
  transformed.

## Value

mrs_data object with additive normally distributed noise.
