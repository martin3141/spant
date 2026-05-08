# Generate trapezoidal regressors.

Generate trapezoidal regressors.

## Usage

``` r
gen_trap_reg(
  onset,
  duration,
  trial_type = NULL,
  mrs_data = NULL,
  tr = NULL,
  Ndyns = NULL,
  Ntrans = NULL,
  rise_t = 0,
  fall_t = 0,
  exp_fall = FALSE,
  exp_fall_power = 1,
  smo_sigma = NULL,
  match_tr = TRUE,
  dt = 0.01,
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

- rise_t:

  time to reach a plateau from baseline in seconds.

- fall_t:

  time to fall from plateau level back to baseline in seconds.

- exp_fall:

  model an exponential fall instead of linear.

- exp_fall_power:

  exponential fall power.

- smo_sigma:

  standard deviation of Gaussian smoothing kernel in seconds. Set to
  NULL to disable (default behavior).

- match_tr:

  match the output to the input mrs_data.

- dt:

  timing resolution for internal calculations.

- normalise:

  normalise the response function to have a maximum value of one.

## Value

trapezoidal regressor data frame.
