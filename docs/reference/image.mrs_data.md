# Image plot method for objects of class mrs_data.

Image plot method for objects of class mrs_data.

## Usage

``` r
# S3 method for class 'mrs_data'
image(
  x,
  xlim = NULL,
  mode = "re",
  col = NULL,
  plot_dim = NULL,
  x_pos = NULL,
  y_pos = NULL,
  z_pos = NULL,
  dyn = 1,
  coil = 1,
  restore_def_par = TRUE,
  y_ticks = NULL,
  hline = NULL,
  hline_lty = 2,
  hline_col = "white",
  vline = NULL,
  vline_lty = 2,
  vline_col = "white",
  legend = FALSE,
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

- col:

  Colour map to use, defaults to viridis.

- plot_dim:

  the dimension to display on the y-axis, can be one of: "dyn",
  "time_sec", x", "y", "z", "coil" or NULL. If NULL (the default) all
  spectra are collapsed into the dynamic dimension and displayed.

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

- restore_def_par:

  restore default plotting par values after the plot has been made.

- y_ticks:

  a vector of indices specifying where to place additional red tick
  marks.

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

- legend:

  add a colour bar to the plot using the imagePlot function from the
  fields package.

- ...:

  other arguments to pass to the plot method.
