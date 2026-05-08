# Read a directory containing Siemens MRS IMA files and combine along the coil dimension. Note that the coil ID is inferred from the sorted file name and should be checked when consistency is required between two directories.

Read a directory containing Siemens MRS IMA files and combine along the
coil dimension. Note that the coil ID is inferred from the sorted file
name and should be checked when consistency is required between two
directories.

## Usage

``` r
read_ima_coil_dir(dir, extra = NULL, verbose = FALSE)
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
