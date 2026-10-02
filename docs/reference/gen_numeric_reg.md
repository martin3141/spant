# Generate a regressor from a numeric vector.

Generate a regressor from a numeric vector.

## Usage

``` r
gen_numeric_reg(
  in_vec,
  name,
  mrs_data = NULL,
  tr = NULL,
  Ndyns = NULL,
  Ntrans = NULL
)
```

## Arguments

- in_vec:

  a numeric input vector the same length as the number of stored dynamic
  scans.

- name:

  a character vector representing the regressor name.

- mrs_data:

  mrs_data object for timing information.

- tr:

  repetition time.

- Ndyns:

  number of dynamic scans stored, potentially less than Ntrans if block
  averaging has been performed.

- Ntrans:

  number of dynamic scans acquired.

## Value

a data frame with a time column and a single regressor column named
according to the name argument.
