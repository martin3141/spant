# Generate polynomial regressors.

Generate polynomial regressors.

## Usage

``` r
gen_poly_reg(degree, mrs_data = NULL, tr = NULL, Ndyns = NULL, Ntrans = NULL)
```

## Arguments

- degree:

  the degree of the polynomial.

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

polynomial regressors.
