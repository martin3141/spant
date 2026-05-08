# Plot a slice from a 7 dimensional array.

Plot a slice from a 7 dimensional array.

## Usage

``` r
plot_slice_map(
  data,
  zlim = NULL,
  mask_map = NULL,
  mask_cutoff = 20,
  interp = 1,
  slice = 1,
  dyn = 1,
  coil = 1,
  ref = 1,
  denom = NULL,
  horizontal = FALSE
)
```

## Arguments

- data:

  7d array of values to be plotted.

- zlim:

  smallest and largest values to be plotted.

- mask_map:

  matching map with logical values to indicate if the corresponding
  values should be plotted.

- mask_cutoff:

  minimum values to plot (as a percentage of the maximum).

- interp:

  map interpolation factor.

- slice:

  the slice index to plot.

- dyn:

  the dynamic index to plot.

- coil:

  the coil element number to plot.

- ref:

  reference index to plot.

- denom:

  map to use as a denominator.

- horizontal:

  display the colourbar horizontally (logical).
