# Generate impulse regressors.

Generate impulse regressors.

## Usage

``` r
gen_impulse_reg(
  onset,
  trial_type = NULL,
  mrs_data = NULL,
  tr = NULL,
  Ndyns = NULL,
  Ntrans = NULL
)
```

## Arguments

- onset:

  stimulus onset in seconds.

- trial_type:

  string label for the stimulus.

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

impulse regressors data frame.
