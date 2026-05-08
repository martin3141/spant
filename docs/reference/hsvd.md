# HSVD of an mrs_data object.

HSVD method as described in: Barkhuijsen H, de Beer R, van Ormondt D.
Improved algorithm for noniterative and timedomain model fitting to
exponentially damped magnetic resonance signals. J Magn Reson
1987;73:553-557.

## Usage

``` r
hsvd(mrs_data, comps = 40, irlba = TRUE, max_damp = 10)
```

## Arguments

- mrs_data:

  mrs_data object to be decomposed.

- comps:

  number of Lorentzian components to use for modelling.

- irlba:

  option to use irlba SVD (logical).

- max_damp:

  maximum allowable damping factor.

## Value

basis matrix and signal table.
