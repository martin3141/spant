# Generate BOLD regressors.

Generate BOLD regressors.

## Usage

``` r
gen_bold_reg(
  onset,
  duration = NULL,
  trial_type = NULL,
  mrs_data = NULL,
  tr = NULL,
  Ndyns = NULL,
  Ntrans = NULL,
  match_tr = TRUE,
  dt = 0.1,
  normalise = FALSE
)
```

## Arguments

- onset:

  stimulus onset in seconds.

- duration:

  stimulus duration in seconds.

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

- match_tr:

  match the output to the input mrs_data.

- dt:

  timing resolution for internal calculations.

- normalise:

  normalise the response function to have a maximum value of one.

## Value

BOLD regressor data frame.
