# Perform a mathematical operation on a spectral region.

Perform a mathematical operation on a spectral region.

## Usage

``` r
spec_op(
  mrs_data,
  xlim = NULL,
  operator = "sum",
  freq_scale = "ppm",
  mode = "re"
)
```

## Arguments

- mrs_data:

  MRS data.

- xlim:

  spectral range to be integrated (defaults to full range).

- operator:

  can be "sum" (default), "mean", "l2", "max", "max_cplx, "min" or
  "max-min".

- freq_scale:

  units of xlim, can be : "ppm", "hz" or "points".

- mode:

  spectral mode, can be : "re", "im", "mod" or "cplx".

## Value

an array of integral values.
