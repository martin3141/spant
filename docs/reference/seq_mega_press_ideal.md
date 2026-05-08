# MEGA-PRESS sequence with ideal localisation pulses and Gaussian shaped editing pulse.

MEGA-PRESS sequence with ideal localisation pulses and Gaussian shaped
editing pulse.

## Usage

``` r
seq_mega_press_ideal(
  spin_params,
  ft,
  ref,
  ed_freq = 1.89,
  TE1 = 0.015,
  TE2 = 0.053,
  BW = 110,
  steps = 50
)
```

## Arguments

- spin_params:

  spin system definition.

- ft:

  transmitter frequency in Hz.

- ref:

  reference value for ppm scale.

- ed_freq:

  editing pulse frequency in ppm.

- TE1:

  TE1 sequence parameter in seconds (TE=TE1+TE2).

- TE2:

  TE2 sequence parameter in seconds.

- BW:

  editing pulse bandwidth in Hz.

- steps:

  number of hard pulses used to approximate the editing pulse.

## Value

list of resonance amplitudes and frequencies.
