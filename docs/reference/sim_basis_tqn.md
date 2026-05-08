# Simulate a basis file using TARQUIN.

Simulate a basis file using TARQUIN.

## Usage

``` r
sim_basis_tqn(
  fs = def_fs(),
  ft = def_ft(),
  N = def_N(),
  ref = def_ref(),
  opts = NULL
)
```

## Arguments

- fs:

  sampling frequency

- ft:

  transmitter frequency

- N:

  number of data points

- ref:

  chemical shift reference

- opts:

  list of options to pass to TARQUIN.

## Examples

``` r
if (FALSE) { # \dontrun{
write_basis_tqn('test.basis',mrs_data,c("--echo","0.04"))
} # }
```
