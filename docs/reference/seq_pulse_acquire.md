# Simple pulse and acquire sequence with ideal pulses.

Simple pulse and acquire sequence with ideal pulses.

## Usage

``` r
seq_pulse_acquire(spin_params, ft, ref, nuc = "1H", acq_delay = 0)
```

## Arguments

- spin_params:

  spin system definition.

- ft:

  transmitter frequency in Hz.

- ref:

  reference value for ppm scale.

- nuc:

  acquisition nucleus.

- acq_delay:

  delay between excitation and acquisition.

## Value

list of resonance amplitudes and frequencies.
