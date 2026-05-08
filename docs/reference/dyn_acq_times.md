# Return a time scale vector of acquisition times for a dynamic MRS scan. The first temporal scan is assigned a value of 0.

Return a time scale vector of acquisition times for a dynamic MRS scan.
The first temporal scan is assigned a value of 0.

## Usage

``` r
dyn_acq_times(mrs_data = NULL, tr = NULL, Ndyns = NULL, Ntrans = NULL)
```

## Arguments

- mrs_data:

  MRS data.

- tr:

  repetition time.

- Ndyns:

  number of dynamic scans stored, potentially less than Ntrans if block
  averaging has been performed.

- Ntrans:

  number of dynamic scans acquired.

## Value

time scale vector in units of seconds.
