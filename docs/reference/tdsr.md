# Time-domain spectral registration.

An implementation of the method published by Near et al MRM 73:44-50
(2015).

## Usage

``` r
tdsr(mrs_data, ref = NULL, xlim = c(4, 0.5), max_t = 0.2)
```

## Arguments

- mrs_data:

  MRS data to be corrected.

- ref:

  optional MRS data to use as a reference, the mean of all dynamics is
  used if this argument is not supplied.

- xlim:

  optional frequency range to perform optimisation, set to NULL to use
  the full range.

- max_t:

  truncate the FID when longer than max_t to reduce time taken.

## Value

a list containing the corrected data; phase and shift values in units of
degrees and Hz respectively.
