# Write MRS data object to file.

Write MRS data object to file.

## Usage

``` r
write_mrs(mrs_data, fname, format = NULL, force = FALSE)
```

## Arguments

- mrs_data:

  object to be written to file, or list of mrs_data objects.

- fname:

  one or more filenames to output.

- format:

  string describing the data format. Must be one of the following :
  "nifti", "dpt", "lcm_raw", "rds". If not specified, the format will be
  guessed from the filename extension.

- force:

  set to TRUE to overwrite any existing files.
