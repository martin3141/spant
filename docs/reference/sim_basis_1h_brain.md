# Simulate a basis-set suitable for 1H brain MRS analysis acquired with a PRESS sequence. Note, ideal pulses are assumed.

Simulate a basis-set suitable for 1H brain MRS analysis acquired with a
PRESS sequence. Note, ideal pulses are assumed.

## Usage

``` r
sim_basis_1h_brain(
  pul_seq = seq_press_ideal,
  acq_paras = def_acq_paras(),
  xlim = c(0.5, 4.2),
  lcm_compat = FALSE,
  ...
)
```

## Arguments

- pul_seq:

  pulse sequence function to use.

- acq_paras:

  list of acquisition parameters or an mrs_data object. See
  [`def_acq_paras`](https://martin3141.github.io/spant/reference/def_acq_paras.md).

- xlim:

  range of frequencies to simulate in ppm.

- lcm_compat:

  exclude lipid and MM signals for use with default LCModel options.

- ...:

  extra parameters to pass to the pulse sequence function.

## Value

basis object.
