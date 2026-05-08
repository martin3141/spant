# PRESS sequence with shaped refocusing pulses.

PRESS sequence with shaped refocusing pulses.

## Usage

``` r
seq_press_2d_shaped(
  spin_params,
  ft,
  ref,
  TE1 = 0.01,
  TE2 = 0.02,
  pulse_file,
  pulse_dur,
  pulse_file_format,
  refoc_flip_angle = 180,
  xy_pulse_ppm = NULL,
  resamp = TRUE,
  fs_resamp = 1e-04
)
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

- pulse_file:

  path to refocusing pulse file.

- pulse_dur:

  refocusing pulse duration.

- pulse_file_format:

  file format for the refocusing pulse.

- refoc_flip_angle:

  refocusing pulse flip angle in degrees (defaults to 180).

- xy_pulse_ppm:

  a vector of ppm values for the offset of each sub-simulation.

- resamp:

  option to resample the pulse shape.

- fs_resamp:

  sampling frequency (Hz) to resample.

## Value

list of resonance amplitudes and frequencies.
