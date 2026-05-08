# Generate regressors by convolving a specified response function with a stimulus.

Generate regressors by convolving a specified response function with a
stimulus.

## Usage

``` r
gen_conv_reg(
  onset,
  duration = NULL,
  trial_type = NULL,
  mrs_data = NULL,
  tr = NULL,
  Ndyns = NULL,
  Ntrans = NULL,
  resp_fn,
  match_tr = TRUE,
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

- resp_fn:

  a data frame specifying the response function to be convolved.

- match_tr:

  match the output to the input mrs_data.

- normalise:

  normalise the response function to have a maximum value of one.

## Value

BOLD regressor data frame.
