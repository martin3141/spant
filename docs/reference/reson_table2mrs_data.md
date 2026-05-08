# Generate mrs_data from a table of single Lorentzian resonances.

Generate mrs_data from a table of single Lorentzian resonances.

## Usage

``` r
reson_table2mrs_data(
  reson_table,
  acq_paras = def_acq_paras(),
  back_extrap_pts = 0
)
```

## Arguments

- reson_table:

  as produced by the hsvd function.

- acq_paras:

  list of acquisition parameters. See

- back_extrap_pts:

  number of data points to back extrapolate
  [`def_acq_paras`](https://martin3141.github.io/spant/reference/def_acq_paras.md)

## Value

mrs_data object.
