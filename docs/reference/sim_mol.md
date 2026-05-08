# Simulate a `mol_parameter` object.

Simulate a `mol_parameter` object.

## Usage

``` r
sim_mol(
  mol,
  pul_seq = seq_pulse_acquire,
  ft = def_ft(),
  ref = def_ref(),
  fs = def_fs(),
  N = def_N(),
  xlim = NULL,
  ...
)
```

## Arguments

- mol:

  `mol_parameter` object.

- pul_seq:

  pulse sequence function to use.

- ft:

  transmitter frequency in Hz.

- ref:

  reference value for ppm scale.

- fs:

  sampling frequency in Hz.

- N:

  number of data points in the spectral dimension.

- xlim:

  ppm range limiting signals to be simulated.

- ...:

  extra parameters to pass to the pulse sequence function.

## Value

`mrs_data` object.
