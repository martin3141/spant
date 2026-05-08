# Crop `basis_set` object based on a frequency range.

Crop `basis_set` object based on a frequency range.

## Usage

``` r
crop_basis(basis, xlim = c(4, 0.2), scale = "ppm")
```

## Arguments

- basis:

  basis_set object to be cropped in the spectral dimension.

- xlim:

  range of values to crop in the spectral dimension eg xlim = c(4, 0.2).

- scale:

  the units to use for the frequency scale, can be one of: "ppm", "hz"
  or "points".

## Value

cropped `mrs_data` object.
