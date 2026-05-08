# Create a two dimensional Gaussian window function stored as a matrix.

Create a two dimensional Gaussian window function stored as a matrix.

## Usage

``` r
gausswin_2d(xN, yN, x0, y0, xw, yw)
```

## Arguments

- xN:

  number of pixels in the x dimension.

- yN:

  number of pixels in the y dimension.

- x0:

  centre of window function in the x direction in units of pixels. Note,
  only integer values are applied.

- y0:

  centre of window function in the y direction in units of pixels. Note,
  only integer values are applied.

- xw:

  the reciprocal of the standard deviation of the Gaussian window in x
  direction.

- yw:

  the reciprocal of the standard deviation of the Gaussian window in y
  direction.

## Value

matrix with dimensions fov_yN x fov_xN.
