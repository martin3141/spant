# Resample a basis-set to match a mrs_data acquisition.

Resample a basis-set to match a mrs_data acquisition.

## Usage

``` r
resample_basis(basis, mrs_data, ref_freq_match = TRUE)
```

## Arguments

- basis:

  the basis to be resampled.

- mrs_data:

  the mrs_data to match the number of data points and sampling
  frequency.

- ref_freq_match:

  apply a frequency shift to the basis to match the reference frequency
  (usually 4.65 or 4.68) of the mrs_data.

## Value

resampled basis set object.
