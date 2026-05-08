# Crop `mrs_data` object based on a frequency range.

Crop `mrs_data` object based on a frequency range.

## Usage

``` r
crop_spec(mrs_data, xlim = c(4, 0.2), scale = "ppm")
```

## Arguments

- mrs_data:

  MRS data.

- xlim:

  range of values to crop in the spectral dimension eg xlim = c(4, 0.2).

- scale:

  the units to use for the frequency scale, can be one of: "ppm", "hz"
  or "points".

## Value

cropped `mrs_data` object.
