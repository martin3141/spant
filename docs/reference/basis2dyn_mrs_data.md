# Convert a basis object to a dynamic mrs_data object.

Convert a basis object to a dynamic mrs_data object.

## Usage

``` r
basis2dyn_mrs_data(basis, amps, tr)
```

## Arguments

- basis:

  basis set object.

- amps:

  a data frame with each column corresponding to a basis element and
  each row corresponding to each dynamic scan.

- tr:

  the dataset repetition time in seconds.

## Value

a dynamic mrs_data object.
