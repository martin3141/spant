# Simulate an ideal pulse excitation profile by smoothing a top-hat function with a Gaussian.

Simulate an ideal pulse excitation profile by smoothing a top-hat
function with a Gaussian.

## Usage

``` r
sim_th_excit_profile(bw = 1500, sigma = 50, fa = 180)
```

## Arguments

- bw:

  top-hat bandwidth (Hz).

- sigma:

  Gaussian width smoothing parameter (Hz).

- fa:

  intended flip angle of the pulse.

## Value

data frame containing the frequency scale, excitation profile and
corresponding flip-angles.
