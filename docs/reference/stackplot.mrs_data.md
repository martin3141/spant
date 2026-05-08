# Stackplot plotting method for objects of class mrs_data.

Stackplot plotting method for objects of class mrs_data.

## Usage

``` r
# S3 method for class 'mrs_data'
stackplot(
  x,
  xlim = NULL,
  mode = "re",
  x_units = NULL,
  fd = TRUE,
  col = NULL,
  alpha = NULL,
  x_offset = 0,
  y_offset = 0,
  plot_dim = NULL,
  x_pos = NULL,
  y_pos = NULL,
  z_pos = NULL,
  dyn = 1,
  coil = 1,
  bty = NULL,
  labels = NULL,
  lab_cex = 1,
  bl_lty = NULL,
  bl_lwd = 0.5,
  restore_def_par = TRUE,
  show_grid = NULL,
  grid_nx = NULL,
  grid_ny = NA,
  lwd = NULL,
  vline = NULL,
  vline_lty = 2,
  vline_col = "red",
  mar = NULL,
  y_scale = FALSE,
  ...
)
```

## Arguments

- x:

  object of class mrs_data.

- xlim:

  the range of values to display on the x-axis, eg xlim = c(4,1).

- mode:

  representation of the complex numbers to be plotted, can be one of:
  "re", "im", "mod" or "arg".

- x_units:

  the units to use for the x-axis, can be one of: "ppm", "hz", "points"
  or "seconds".

- fd:

  display data in the frequency-domain (default), or time-domain
  (logical).

- col:

  set the colour of the line, eg col = rgb(1, 0, 0, 0.5).

- alpha:

  set the line transparency, eg alpha = 0.5 is 50% transparency.
  Overrides any transparency levels set by col.

- x_offset:

  separate plots in the x-axis direction by this value. Default value is
  0.

- y_offset:

  separate plots in the y-axis direction by this value.

- plot_dim:

  the dimension to display on the y-axis, can be one of: "dyn", "x",
  "y", "z", "coil" or NULL. If NULL (the default) all spectra will be
  collapsed into the dynamic dimension and displayed.

- x_pos:

  the x index to plot.

- y_pos:

  the y index to plot.

- z_pos:

  the z index to plot.

- dyn:

  the dynamic index to plot.

- coil:

  the coil element number to plot.

- bty:

  option to draw a box around the plot. See ?par.

- labels:

  add labels to each data item.

- lab_cex:

  label size.

- bl_lty:

  linetype for the y = 0 baseline trace. A default value NULL results in
  no baseline being plotted.

- bl_lwd:

  linewith for the y = 0 baseline trace. Defaults to 0.5.

- restore_def_par:

  restore default plotting par values after the plot has been made.

- show_grid:

  plot gridlines behind the data (logical). Defaults to TRUE.

- grid_nx:

  number of cells of the grid in x and y direction. When NULL the grid
  aligns with the tick marks on the corresponding default axis (i.e.,
  tickmarks as computed by axTicks). When NA, no grid lines are drawn in
  the corresponding direction.

- grid_ny:

  as above.

- lwd:

  plot linewidth.

- vline:

  x-value to draw a vertical line.

- vline_lty:

  linetype for the vertical line.

- vline_col:

  colour for the vertical line.

- mar:

  option to adjust the plot margins. See ?par.

- y_scale:

  option to display the y-axis values (logical).

- ...:

  other arguments to pass to the matplot method.
