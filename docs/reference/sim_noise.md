# Simulate an mrs_data object containing simulated Gaussian noise.

Simulate an mrs_data object containing simulated Gaussian noise.

## Usage

``` r
sim_noise(
  sd = 0.1,
  fs = def_fs(),
  ft = def_ft(),
  N = def_N(),
  ref = def_ref(),
  nuc = def_nuc(),
  dyns = 1,
  fd = TRUE
)
```

## Arguments

- sd:

  standard deviation of the noise.

- fs:

  sampling frequency in Hz.

- ft:

  transmitter frequency in Hz.

- N:

  number of data points in the spectral dimension.

- ref:

  reference value for ppm scale.

- nuc:

  resonant nucleus.

- dyns:

  number of dynamic scans to generate.

- fd:

  return data in the frequency-domain (TRUE) or time-domain (FALSE)

## Value

mrs_data object.
