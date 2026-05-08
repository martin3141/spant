# Convert a vector into a mrs_data object.

Convert a vector into a mrs_data object.

## Usage

``` r
vec2mrs_data(
  vec,
  mrs_data = NULL,
  fs = NULL,
  ft = NULL,
  ref = NULL,
  nuc = NULL,
  dyns = 1,
  fd = FALSE
)
```

## Arguments

- vec:

  the data vector.

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

- dyns:

  replicate the data across the dynamic dimension.

- fd:

  flag to indicate if the vector is in the frequency domain (logical).

## Value

mrs_data object.
