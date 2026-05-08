# PRESS sequence with ideal pulses.

PRESS sequence with ideal pulses.

## Usage

``` r
seq_press_ideal(spin_params, ft, ref, TE1 = 0.01, TE2 = 0.02)
```

## Arguments

- spin_params:

  spin system definition.

- ft:

  transmitter frequency in Hz.

- ref:

  reference value for ppm scale.

- TE1:

  TE1 sequence parameter in seconds (TE=TE1+TE2).

- TE2:

  TE2 sequence parameter in seconds.

## Value

list of resonance amplitudes and frequencies.
