# Simulate pulse sequence acquisition.

Simulate pulse sequence acquisition.

## Usage

``` r
acquire(sys, rec_phase = 0, tol = 1e-04, detect = NULL, amp_scale = 1)
```

## Arguments

- sys:

  spin system object.

- rec_phase:

  receiver phase in degrees.

- tol:

  ignore resonance amplitudes below this threshold.

- detect:

  detection nuclei.

- amp_scale:

  scaling factor for the output amplitudes.

## Value

a list of resonance amplitudes and frequencies.
