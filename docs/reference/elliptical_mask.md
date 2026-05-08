# Create an elliptical mask stored as a matrix of logical values.

Create an elliptical mask stored as a matrix of logical values.

## Usage

``` r
elliptical_mask(xN, yN, x0, y0, xr, yr, angle)
```

## Arguments

- xN:

  number of pixels in the x dimension.

- yN:

  number of pixels in the y dimension.

- x0:

  centre of ellipse in the x direction in units of pixels.

- y0:

  centre of ellipse in the y direction in units of pixels.

- xr:

  radius in the x direction in units of pixels.

- yr:

  radius in the y direction in units of pixels.

- angle:

  angle of rotation in degrees.

## Value

logical mask matrix with dimensions fov_yN x fov_xN.
