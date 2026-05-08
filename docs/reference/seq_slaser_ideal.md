# sLASER sequence with ideal pulses.

sLASER sequence with ideal pulses.

## Usage

``` r
seq_slaser_ideal(spin_params, ft, ref, TE1 = 0.008, TE2 = 0.011, TE3 = 0.009)
```

## Arguments

- spin_params:

  spin system definition.

- ft:

  transmitter frequency in Hz.

- ref:

  reference value for ppm scale.

- TE1:

  first echo time (between exc. and 1st echo) in seconds.

- TE2:

  second echo time (between 2nd echo and 4th echo) in seconds.

- TE3:

  third echo time (between 4th echo and 5th echo) in seconds.

## Value

list of resonance amplitudes and frequencies.
