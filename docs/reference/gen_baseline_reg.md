# Generate baseline regressor.

Generate baseline regressor.

## Usage

``` r
gen_baseline_reg(mrs_data = NULL, tr = NULL, Ndyns = NULL, Ntrans = NULL)
```

## Arguments

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

a single baseline regressor with value of 1.
