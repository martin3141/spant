# Simulate an mrs_data object containing complex zero valued samples.

Simulate an mrs_data object containing complex zero valued samples.

## Usage

``` r
sim_zero(
  fs = def_fs(),
  ft = def_ft(),
  N = def_N(),
  ref = def_ref(),
  nuc = def_nuc(),
  dyns = 1
)
```

## Arguments

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

## Value

mrs_data object.
