# Arrange spectral plots in a grid.

Arrange spectral plots in a grid.

## Usage

``` r
# S3 method for class 'mrs_data'
gridplot(
  x,
  rows = NA,
  cols = NA,
  mar = c(0, 0, 0, 0),
  oma = c(3.5, 1, 1, 1),
  bty = "o",
  restore_def_par = TRUE,
  ...
)
```

## Arguments

- x:

  object of class mrs_data.

- rows:

  number of grid rows.

- cols:

  number of grid columns.

- mar:

  option to adjust the plot margins. See ?par.

- oma:

  outer margin area.

- bty:

  option to draw a box around the plot. See ?par.

- restore_def_par:

  restore default plotting par values after the plot has been made.

- ...:

  other arguments to pass to the plot method.
