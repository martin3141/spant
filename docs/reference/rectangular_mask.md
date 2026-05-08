# Create a rectangular mask stored as a matrix of logical values.

Create a rectangular mask stored as a matrix of logical values.

## Usage

``` r
rectangular_mask(xN, yN, x0, y0, xw, yw, angle)
```

## Arguments

- xN:

  number of pixels in the x dimension.

- yN:

  number of pixels in the y dimension.

- x0:

  centre of rectangle in the x direction in units of pixels.

- y0:

  centre of rectangle in the y direction in units of pixels.

- xw:

  width in the x direction in units of pixels.

- yw:

  width in the y direction in units of pixels.

- angle:

  angle of rotation in degrees.

## Value

logical mask matrix with dimensions fov_yN x fov_xN.
