# Extract a subset of MRS data.

Extract a subset of MRS data.

## Usage

``` r
get_subset(
  mrs_data,
  x_set = NULL,
  y_set = NULL,
  z_set = NULL,
  dyn_set = NULL,
  coil_set = NULL,
  fd_set = NULL,
  td_set = NULL
)
```

## Arguments

- mrs_data:

  MRS data object.

- x_set:

  x indices to include in the output (default all).

- y_set:

  y indices to include in the output (default all).

- z_set:

  z indices to include in the output (default all).

- dyn_set:

  dynamic indices to include in the output (default all).

- coil_set:

  coil indices to include in the output (default all).

- fd_set:

  frequency domain data indices to include in the output (default all).

- td_set:

  time-domain indices to include in the output (default all).

## Value

selected subset of MRS data.
