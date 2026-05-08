# Zero-fill MRSI data in the k-space x-y direction.

Zero-fill MRSI data in the k-space x-y direction.

## Usage

``` r
zf_xy(mrs_data, factor = 2)
```

## Arguments

- mrs_data:

  MRSI data.

- factor:

  zero-filling factor, a factor of 2 returns a dataset with twice the
  original points in the x-y directions. Factors smaller than one are
  permitted, such that a factor of 0.5 returns half the k-space points
  in the x-y directions.

## Value

zero-filled data.
