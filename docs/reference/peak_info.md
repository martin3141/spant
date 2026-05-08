# Search for the highest peak in a spectral region and return the frequency, height and FWHM.

Search for the highest peak in a spectral region and return the
frequency, height and FWHM.

## Usage

``` r
peak_info(
  mrs_data,
  xlim = c(4, 0.5),
  interp_f = 4,
  scale = "ppm",
  mode = "real"
)
```

## Arguments

- mrs_data:

  an object of class `mrs_data`.

- xlim:

  frequency range (default units of PPM) to search for the highest peak.

- interp_f:

  interpolation factor, defaults to 4x.

- scale:

  the units to use for the frequency scale, can be one of: "ppm", "hz"
  or "points".

- mode:

  spectral mode, can be : "real", "imag" or "mod".

## Value

list of arrays containing the highest peak frequency, height and FWHM in
units of PPM and Hz.
