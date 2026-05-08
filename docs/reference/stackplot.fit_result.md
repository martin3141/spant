# Plot the fitting results of an object of class `fit_result` with individual basis set components shown.

Plot the fitting results of an object of class `fit_result` with
individual basis set components shown.

## Usage

``` r
# S3 method for class 'fit_result'
stackplot(
  x,
  xlim = NULL,
  y_offset = 0,
  dyn = 1,
  x_pos = 1,
  y_pos = 1,
  z_pos = 1,
  coil = 1,
  n = NULL,
  sub_bl = FALSE,
  labels = FALSE,
  label_names = NULL,
  sig_col = "black",
  restore_def_par = TRUE,
  omit_signals = NULL,
  combine_lipmm = FALSE,
  combine_metab = FALSE,
  mar = NULL,
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

- xlim:

  the range of values to display on the x-axis, eg xlim = c(4,1).

- y_offset:

  separate basis signals in the y-axis direction by this value.

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

- n:

  single index element to plot (overrides other indices when given).

- sub_bl:

  subtract the baseline from the data and fit (logical).

- labels:

  print signal labels at the right side of the plot.

- label_names:

  provide a character vector of signal names to replace the defaults
  determined from the basis set.

- sig_col:

  colour of individual signal components.

- restore_def_par:

  restore default plotting par values after the plot has been made.

- omit_signals:

  a character vector of basis signal names to be removed from the plot.

- combine_lipmm:

  combine all basis signals with names starting with "Lip" or "MM".

- combine_metab:

  combine all basis signals with names not starting with "Lip" or "MM".

- mar:

  option to adjust the plot margins. See ?par.

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
