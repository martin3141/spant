# Plot the fitting results of an object of class `fit_result`.

Plot the fitting results of an object of class `fit_result`.

## Usage

``` r
# S3 method for class 'fit_result'
plot(
  x,
  dyn = 1,
  x_pos = 1,
  y_pos = 1,
  z_pos = 1,
  coil = 1,
  xlim = NULL,
  data_only = FALSE,
  label = NULL,
  plot_sigs = NULL,
  n = NULL,
  sub_bl = FALSE,
  mar = NULL,
  restore_def_par = TRUE,
  ylim = NULL,
  y_scale = FALSE,
  show_grid = TRUE,
  grid_nx = NULL,
  grid_ny = NA,
  invert_fit = FALSE,
  ...
)
```

## Arguments

- x:

  fit_result object.

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

- xlim:

  the range of values to display on the x-axis, eg xlim = c(4,1).

- data_only:

  display only the processed data (logical).

- label:

  character string to add to the top left of the plot window.

- plot_sigs:

  a character vector of signal names to add to the plot.

- n:

  single index element to plot (overrides other indices when given).

- sub_bl:

  subtract the baseline from the data and fit (logical).

- mar:

  option to adjust the plot margins. See ?par.

- restore_def_par:

  restore default plotting par values after the plot has been made.

- ylim:

  range of values to display on the y-axis, eg ylim = c(0,10).

- y_scale:

  option to display the y-axis values (logical).

- show_grid:

  plot gridlines behind the data (logical). Defaults to TRUE.

- grid_nx:

  number of cells of the grid in x and y direction. When NULL the grid
  aligns with the tick marks on the corresponding default axis (i.e.,
  tickmarks as computed by axTicks). When NA, no grid lines are drawn in
  the corresponding direction.

- grid_ny:

  as above.

- invert_fit:

  show the fit result "upside-down"/

- ...:

  further arguments to plot method.
