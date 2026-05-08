# Extract the fit amplitudes from an object of class `fit_result`.

Extract the fit amplitudes from an object of class `fit_result`.

## Usage

``` r
fit_amps(
  x,
  inc_index = FALSE,
  sort_names = FALSE,
  append_common_1h_comb = TRUE
)
```

## Arguments

- x:

  `fit_result` object.

- inc_index:

  include columns for the voxel index.

- sort_names:

  sort the basis set names alphabetically.

- append_common_1h_comb:

  append commonly used 1H metabolite combinations eg tNAA = NAA + NAAG.

## Value

a dataframe of amplitudes.
