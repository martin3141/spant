# Scale mrs_data to a spectral region.

Scale mrs_data to a spectral region.

## Usage

``` r
scale_spec(
  mrs_data,
  xlim = NULL,
  operator = "sum",
  freq_scale = "ppm",
  mode = "re",
  mean_dyns = NULL,
  ret_scale_factor = FALSE
)
```

## Arguments

- mrs_data:

  MRS data.

- xlim:

  spectral range to be integrated (defaults to full range).

- operator:

  can be "sum" (default), "mean", "l2", "max", "min" or "max-min".

- freq_scale:

  units of xlim, can be : "ppm", "Hz" or "points".

- mode:

  spectral mode, can be : "re", "im", "mod" or "cplx".

- mean_dyns:

  mean the dynamic scans before applying the operator. The same scaling
  value will be applied to each individual dynamic.

- ret_scale_factor:

  option to return the scaling factor in addition to the scaled data.

## Value

normalised data.
