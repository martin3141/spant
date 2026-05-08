# Plotting method for objects of class mrs_data.

Plotting method for objects of class mrs_data.

## Usage

``` r
# S3 method for class 'mrs_data'
plot(
  x,
  dyn = 1,
  x_pos = 1,
  y_pos = 1,
  z_pos = 1,
  coil = 1,
  fd = TRUE,
  x_units = NULL,
  xlim = NULL,
  y_scale = FALSE,
  x_ax = TRUE,
  mode = "re",
  lwd = NULL,
  bty = NULL,
  label = "",
  restore_def_par = TRUE,
  mar = NULL,
  xaxis_lab = NULL,
  yaxis_lab = NULL,
  xat = NULL,
  xlabs = TRUE,
  yat = NULL,
  ylabs = TRUE,
  show_grid = TRUE,
  grid_nx = NULL,
  grid_ny = NA,
  col = NULL,
  alpha = NULL,
  bl_lty = NULL,
  hline = NULL,
  hline_lty = 2,
  hline_col = "red",
  vline = NULL,
  vline_lty = 2,
  vline_col = "red",
  ...
)
```

## Arguments

- x:

  object of class mrs_data.

- dyn:

  the dynamic index to plot.

- x_pos:

  the x index to plot.

- y_pos:

  the y index to plot.

- z_pos:

  the z index to plot.

- coil:

  the coil element number to plot.

- fd:

  display data in the frequency-domain (default), or time-domain
  (logical).

- x_units:

  the units to use for the x-axis, can be one of: "ppm", "hz", "points",
  "seconds" or "ms".

- xlim:

  the range of values to display on the x-axis, eg xlim = c(4,1).

- y_scale:

  option to display the y-axis values (logical).

- x_ax:

  option to display the x-axis values (logical).

- mode:

  representation of the complex numbers to be plotted, can be one of:
  "re", "im", "mod" or "arg".

- lwd:

  plot linewidth.

- bty:

  option to draw a box around the plot. See ?par.

- label:

  character string to add to the top left of the plot window.

- restore_def_par:

  restore default plotting par values after the plot has been made.

- mar:

  option to adjust the plot margins. See ?par.

- xaxis_lab:

  x-axis label.

- yaxis_lab:

  y-axis label.

- xat:

  x-axis tick label values.

- xlabs:

  x-axis tick labels.

- yat:

  y-axis tick label values.

- ylabs:

  y-axis tick labels.

- show_grid:

  plot gridlines behind the data (logical). Defaults to TRUE.

- grid_nx:

  number of cells of the grid in x and y direction. When NULL the grid
  aligns with the tick marks on the corresponding default axis (i.e.,
  tickmarks as computed by axTicks). When NA, no grid lines are drawn in
  the corresponding direction.

- grid_ny:

  as above.

- col:

  set the line colour, eg col = rgb(0.5, 0.5, 0.5).

- alpha:

  set the line transparency, eg alpha = 0.5 is 50% transparency.
  Overrides any transparency levels set by col.

- bl_lty:

  linetype for the y = 0 baseline trace. A default value NULL results in
  no baseline being plotted.

- hline:

  add a horizontal line at the specified value.

- hline_lty:

  linetype for the horizontal line.

- hline_col:

  colour for the horizontal line.

- vline:

  add a vertical line at the specified value.

- vline_lty:

  linetype for the vertical line.

- vline_col:

  colour for the vertical line.

- ...:

  other arguments to pass to the plot method.
