# Convert default LCM/TARQUIN concentration scaling to molal units with partial volume correction.

Convert default LCM/TARQUIN concentration scaling to molal units with
partial volume correction.

## Usage

``` r
scale_amp_molar2molal_pvc(fit_result, p_vols, te, tr)
```

## Arguments

- fit_result:

  a `fit_result` object to apply partial volume correction.

- p_vols:

  a numeric vector of partial volumes expressed as percentages. For
  example, a voxel containing 100% white matter tissue would use :
  p_vols = c(WM = 100, GM = 0, CSF = 0).

- te:

  the MRS TE.

- tr:

  the MRS TR.

## Value

a `fit_result` object with a rescaled results table.
