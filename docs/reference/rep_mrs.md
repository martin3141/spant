# Replicate a scan over a given dimension.

Replicate a scan over a given dimension.

## Usage

``` r
rep_mrs(
  mrs_data,
  x_rep = 1,
  y_rep = 1,
  z_rep = 1,
  dyn_rep = 1,
  coil_rep = 1,
  warn = TRUE
)
```

## Arguments

- mrs_data:

  MRS data to be replicated.

- x_rep:

  number of x replications.

- y_rep:

  number of y replications.

- z_rep:

  number of z replications.

- dyn_rep:

  number of dynamic replications.

- coil_rep:

  number of coil replications.

- warn:

  print a warning when the data dimensions do not change.

## Value

replicated data object.
