# Integrate a spectral region.

See spec_op function for a more complete set of spectral operations.

## Usage

``` r
int_spec(mrs_data, xlim = NULL, freq_scale = "ppm", mode = "re")
```

## Arguments

- mrs_data:

  MRS data.

- xlim:

  spectral range to be integrated (defaults to full range).

- freq_scale:

  units of xlim, can be : "ppm", "hz" or "points".

- mode:

  spectral mode, can be : "re", "im", "mod" or "cplx".

## Value

an array of integral values.
