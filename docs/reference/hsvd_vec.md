# HSVD of a complex vector.

HSVD method as described in: Barkhuijsen H, de Beer R, van Ormondt D.
Improved algorithm for noniterative and timedomain model fitting to
exponentially damped magnetic resonance signals. J Magn Reson
1987;73:553-557.

## Usage

``` r
hsvd_vec(y, fs, comps = 40, irlba = TRUE, max_damp = 0)
```

## Arguments

- y:

  time domain signal to be filtered as a vector.

- fs:

  sampling frequency of y.

- comps:

  number of Lorentzian components to use for modelling.

- irlba:

  option to use irlba SVD (logical).

- max_damp:

  maximum allowable damping factor. Default value of 0 ensures resultant
  model is damped.

## Value

basis matrix and signal table.
