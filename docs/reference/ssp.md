# Signal space projection method for lipid suppression.

Signal space projection method as described in: Tsai SY, Lin YR, Lin HY,
Lin FH. Reduction of lipid contamination in MR spectroscopy imaging
using signal space projection. Magn Reson Med 2019 Mar;81(3):1486-1498.

## Usage

``` r
ssp(mrs_data, comps = 5, xlim = c(1.5, 0.8))
```

## Arguments

- mrs_data:

  MRS data object.

- comps:

  the number of spatial components to use.

- xlim:

  spectral range (in ppm) covering the lipid signals.

## Value

lipid suppressed `mrs_data` object.
