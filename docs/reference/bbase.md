# Generate a spline basis, slightly adapted from : "Splines, knots, and penalties", Eilers 2010.

Generate a spline basis, slightly adapted from : "Splines, knots, and
penalties", Eilers 2010.

## Usage

``` r
bbase(N, number, deg = 3)
```

## Arguments

- N:

  number of data points.

- number:

  number of spline functions.

- deg:

  spline degree : deg = 1 linear, deg = 2 quadratic, deg = 3 cubic.

## Value

spline basis as a matrix.
