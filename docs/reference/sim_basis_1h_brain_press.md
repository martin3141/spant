# Simulate a basis-set suitable for 1H brain MRS analysis acquired with a PRESS sequence. Note, ideal pulses are assumed.

Simulate a basis-set suitable for 1H brain MRS analysis acquired with a
PRESS sequence. Note, ideal pulses are assumed.

## Usage

``` r
sim_basis_1h_brain_press(
  acq_paras = def_acq_paras(),
  xlim = c(0.5, 4.2),
  lcm_compat = FALSE,
  TE1 = 0.01,
  TE2 = 0.02
)
```

## Arguments

- acq_paras:

  list of acquisition parameters or an mrs_data object. See
  [`def_acq_paras`](https://martin3141.github.io/spant/reference/def_acq_paras.md)

- xlim:

  range of frequencies to simulate in ppm.

- lcm_compat:

  exclude lipid and MM signals for use with default LCModel options.

- TE1:

  TE1 of PRESS sequence (TE = TE1 + TE2).

- TE2:

  TE2 of PRESS sequence.

## Value

basis object.
