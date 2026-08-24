# Generate spline regressors.

Generate spline regressors.

## Usage

``` r
gen_spline_reg(number, mrs_data = NULL, tr = NULL, Ndyns = NULL, Ntrans = NULL)
```

## Arguments

- number:

  the number of spline functions to generate.

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

spline regressors.
