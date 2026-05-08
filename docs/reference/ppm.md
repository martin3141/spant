# Return the ppm scale of an MRS dataset or fit result.

Return the ppm scale of an MRS dataset or fit result.

## Usage

``` r
ppm(x, ft = NULL, ref = NULL, fs = NULL, N = NULL)

# S3 method for class 'mrs_data'
ppm(x, ft = NULL, ref = NULL, fs = NULL, N = NULL)

# S3 method for class 'fit_result'
ppm(x, ft = NULL, ref = NULL, fs = NULL, N = NULL)
```

## Arguments

- x:

  MRS dataset of fit result.

- ft:

  transmitter frequency in Hz, does not apply when the object is a fit
  result.

- ref:

  reference value for ppm scale, does not apply when the object is a fit
  result.

- fs:

  sampling frequency in Hz, does not apply when the object is a fit
  result.

- N:

  number of data points in the spectral dimension, does not apply when
  the object is a fit result.

## Value

ppm scale.
