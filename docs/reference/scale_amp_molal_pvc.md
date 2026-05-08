# Apply water reference scaling to a fitting results object to yield metabolite quantities in millimolar (mM) units (mol / kg of tissue water).

Details of this method can be found in "Use of tissue water as a
concentration reference for proton spectroscopic imaging" by Gasparovic
et al MRM 2006 55(6):1219-26. 1.5 Tesla relaxation assumptions are taken
from this paper. For 3 Tesla data, relaxation assumptions are taken from
"NMR relaxation times in the human brain at 3.0 Tesla" by Wansapura et
al J Magn Reson Imaging 1999 9(4):531-8.

## Usage

``` r
scale_amp_molal_pvc(fit_result, ref_data, p_vols, te, tr, ...)
```

## Arguments

- fit_result:

  result object generated from fitting.

- ref_data:

  water reference MRS data object.

- p_vols:

  a numeric vector of partial volumes expressed as percentages. For
  example, a voxel containing 100% white matter tissue would use :
  p_vols = c(WM = 100, GM = 0, CSF = 0).

- te:

  the MRS TE in seconds.

- tr:

  the MRS TR in seconds.

- ...:

  additional arguments to get_td_amp function.

## Value

A `fit_result` object with a rescaled results table.
