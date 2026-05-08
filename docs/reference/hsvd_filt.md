# HSVD based signal filter.

HSVD based signal filter described in: Barkhuijsen H, de Beer R, van
Ormondt D. Improved algorithm for noniterative and timedomain model
fitting to exponentially damped magnetic resonance signals. J Magn Reson
1987;73:553-557.

## Usage

``` r
hsvd_filt(
  mrs_data,
  xlim = c(-30, 30),
  comps = 40,
  irlba = TRUE,
  max_damp = 10,
  scale = "hz",
  return_model = FALSE
)
```

## Arguments

- mrs_data:

  MRS data to be filtered.

- xlim:

  frequency range to filter, default units are Hz which can be changed
  to ppm using the "scale" argument.

- comps:

  number of Lorentzian components to use for modelling.

- irlba:

  option to use irlba SVD (logical).

- max_damp:

  maximum allowable damping factor.

- scale:

  either "hz" or "ppm" to set the frequency units of xlim.

- return_model:

  by default the filtered spectrum is returned. Set return_model to TRUE
  to return the HSVD model of the data.

## Value

filtered data or model depending on the return_model argument.
