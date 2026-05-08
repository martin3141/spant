# Convert a 7 dimensional array in into a mrs_data object. The array dimensions should be ordered as : dummy, X, Y, Z, dynamic, coil, FID.

Convert a 7 dimensional array in into a mrs_data object. The array
dimensions should be ordered as : dummy, X, Y, Z, dynamic, coil, FID.

## Usage

``` r
array2mrs_data(
  data_array,
  mrs_data = NULL,
  fs = NULL,
  ft = NULL,
  ref = NULL,
  nuc = NULL,
  fd = FALSE
)
```

## Arguments

- data_array:

  7d data array.

- mrs_data:

  example data to copy acquisition parameters from.

- fs:

  sampling frequency in Hz.

- ft:

  transmitter frequency in Hz.

- ref:

  reference value for ppm scale.

- nuc:

  nucleus that is resonant at the transmitter frequency.

- fd:

  flag to indicate if the matrix is in the frequency domain (logical).

## Value

mrs_data object.
