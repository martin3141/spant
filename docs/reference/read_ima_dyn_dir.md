# Read a directory containing Siemens MRS IMA files and combine along the dynamic dimension. Note that the coil ID is inferred from the sorted file name and should be checked when consistency is required.

Read a directory containing Siemens MRS IMA files and combine along the
dynamic dimension. Note that the coil ID is inferred from the sorted
file name and should be checked when consistency is required.

## Usage

``` r
read_ima_dyn_dir(dir, extra = NULL, verbose = FALSE)
```

## Arguments

- dir:

  data directory path.

- extra:

  an optional data frame to provide additional variables for use in
  subsequent analysis steps, eg id or grouping variables.

- verbose:

  output extra information to the console.

## Value

mrs_data object.
