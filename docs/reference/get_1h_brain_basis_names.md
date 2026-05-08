# Return a character vector of common 1H molecules found in healthy human brain.

Note, this is a basic set and it may be appropriate to also include Asc,
Gly and PEth for high quality MRS data.

## Usage

``` r
get_1h_brain_basis_names(add = NULL, remove = NULL, inc_lip_mm = TRUE)
```

## Arguments

- add:

  optional character vector of additional molecular names. Eg c("asc",
  "gly", "peth").

- remove:

  optional character vector of molecular names to remove from the set.
  Eg c("m_cr_ch2").

- inc_lip_mm:

  include Lipid and MM basis signals.

## Value

a character vector of molecule names.
