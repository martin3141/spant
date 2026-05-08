# STEAM sequence with ideal pulses.

STEAM sequence with ideal pulses.

## Usage

``` r
seq_steam_ideal(spin_params, ft, ref, TE = 0.03, TM = 0.02, amp_scale = 2)
```

## Arguments

- spin_params:

  spin system definition.

- ft:

  transmitter frequency in Hz.

- ref:

  reference value for ppm scale.

- TE:

  sequence parameter in seconds.

- TM:

  sequence parameter in seconds.

- amp_scale:

  amplitude scaling factor. Set to 2 (default) to ensure correct scaling
  for water reference scaling. Set to 1 to maintain the inherent loss of
  signal associated with STEAM.

## Value

list of resonance amplitudes and frequencies.
