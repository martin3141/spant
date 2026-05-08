# Convert a matrix (with spectral points in the column dimension and dynamics in the row dimensions) into a mrs_data object.

Convert a matrix (with spectral points in the column dimension and
dynamics in the row dimensions) into a mrs_data object.

## Usage

``` r
mat2mrs_data(
  mat,
  mrs_data = NULL,
  fs = NULL,
  ft = NULL,
  ref = NULL,
  nuc = NULL,
  fd = FALSE
)
```

## Arguments

- mat:

  data matrix.

- mrs_data:

  example data to copy acquisition parameters from.

- fs:

  sampling frequency in Hz.

- ft:

  transmitter frequency in Hz.

- ref:

  reference value for ppm scale.

- nuc:

  resonant nucleus.

- fd:

  flag to indicate if the matrix is in the frequency domain (logical).

## Value

mrs_data object.
