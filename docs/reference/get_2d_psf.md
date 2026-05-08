# Get the point spread function (PSF) for a 2D phase encoded MRSI scan.

Get the point spread function (PSF) for a 2D phase encoded MRSI scan.

## Usage

``` r
get_2d_psf(
  FOV = 160,
  mat_size = 16,
  sampling = "circ",
  hamming = FALSE,
  ensure_odd = TRUE
)
```

## Arguments

- FOV:

  field of view in mm.

- mat_size:

  acquisition matrix size (not interpolated).

- sampling:

  can be either "circ" for circular or "rect" for rectangular.

- hamming:

  should Hamming k-space weighting be applied (default FALSE).

- ensure_odd:

  add 1mm to the FOV when required to ensure the output pdf has odd
  dimensions. Required when using get_mrsi2d_seg.

## Value

A matrix of the PSF with 1mm resolution.
