# Return a list of `mol_parameter` objects suitable for 1H brain MRS analyses.

Return a list of `mol_parameter` objects suitable for 1H brain MRS
analyses.

## Usage

``` r
get_1h_brain_basis_paras_v3(ft, metab_lw = NULL, lcm_compat = FALSE)
```

## Arguments

- ft:

  transmitter frequency in Hz.

- metab_lw:

  linewidth of metabolite signals (Hz).

- lcm_compat:

  when TRUE, lipid, MM and -CrCH molecules will be excluded from the
  output.

## Value

list of `mol_parameter` objects.
