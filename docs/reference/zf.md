# Zero-fill MRS data in the time domain.

Zero-fill MRS data in the time domain.

## Usage

``` r
zf(x, factor = 2, offset = 0)

# S3 method for class 'list'
zf(x, factor = 2, offset = 0)

# S3 method for class 'mrs_data'
zf(x, factor = 2, offset = 0)

# S3 method for class 'basis_set'
zf(x, factor = 2, offset = 0)
```

## Arguments

- x:

  input mrs_data or basis_set object.

- factor:

  zero-filling factor, factor of 2 returns a dataset with twice the
  original data points.

- offset:

  number of points from the end of the FID to insert the zero values.

## Value

zero-filled data.
