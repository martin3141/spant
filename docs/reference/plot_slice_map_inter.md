# Plot an interactive slice map from a data array where voxels can be selected to display a corresponding spectrum.

Plot an interactive slice map from a data array where voxels can be
selected to display a corresponding spectrum.

## Usage

``` r
plot_slice_map_inter(
  mrs_data,
  map = NULL,
  xlim = NULL,
  slice = 1,
  zlim = NULL,
  mask_map = NULL,
  denom = NULL,
  mask_cutoff = 20,
  interp = 1,
  mode = "re",
  y_scale = FALSE,
  ylim = NULL,
  coil = 1,
  fd = TRUE
)
```

## Arguments

- mrs_data:

  spectral data.

- map:

  array of values to be plotted, defaults to the integration of the
  modulus of the full spectral width.

- xlim:

  spectral region to plot.

- slice:

  the slice index to plot.

- zlim:

  smallest and largest values to be plotted.

- mask_map:

  matching map with logical values to indicate if the corresponding
  values should be plotted.

- denom:

  map to use as a denominator.

- mask_cutoff:

  minimum values to plot (as a percentage of the maximum).

- interp:

  map interpolation factor.

- mode:

  representation of the complex spectrum to be plotted, can be one of:
  "re", "im", "mod" or "arg".

- y_scale:

  option to display the y-axis values (logical).

- ylim:

  intensity range to plot.

- coil:

  coil element to plot.

- fd:

  display data in the frequency-domain (default), or time-domain
  (logical).
