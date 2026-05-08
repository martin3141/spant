# Plot a 2D slice from an MRSI fit result object.

Plot a 2D slice from an MRSI fit result object.

## Usage

``` r
plot_slice_fit(
  fit_res,
  map,
  map_denom = NULL,
  slice = 1,
  zlim = NULL,
  interp = 1
)
```

## Arguments

- fit_res:

  `fit_result` object.

- map:

  fit result values to display as a colour map. Can be specified as a
  character string or array of numeric values. Defaults to "tNAA".

- map_denom:

  fit result values to divide the map argument by. Can be specified as a
  character string (eg "tCr") or array of numeric values.

- slice:

  slice to plot in the z direction.

- zlim:

  range of values to plot.

- interp:

  interpolation factor.
