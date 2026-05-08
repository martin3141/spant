# Apply line-broadening to an mrs_data object to achieve a specified linewidth.

Apply line-broadening to an mrs_data object to achieve a specified
linewidth.

## Usage

``` r
set_lw(mrs_data, lw, xlim = c(4, 0.5), lg = 1, mask_narrow = TRUE)
```

## Arguments

- mrs_data:

  data in.

- lw:

  target linewidth in units of ppm.

- xlim:

  region to search for peaks to obtain a linewidth estimate.

- lg:

  Lorentz-Gauss lineshape parameter.

- mask_narrow:

  masks spectra where the requested linewidth is too narrow, if set
  FALSE the spectra are not changed.

## Value

line-broadened data.
