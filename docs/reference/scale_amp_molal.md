# Apply water reference scaling to a fitting results object to yield metabolite quantities in millimolar (mM) units (mol / kg of tissue water).

Note, this function assumes the volume contains a homogeneous voxel, eg
pure WM, GM or CSF. Also note that in the case of a homogeneous voxel
the relative densities of MR-visible water (eg GM=0.78, WM=0.65, and
CSF=0.97) cancel out and don't need to be considered. Use
scale_amp_molal_pvc for volumes containing multiple compartments.
Details of this method can be found in "Use of tissue water as a
concentration reference for proton spectroscopic imaging" by Gasparovic
et al MRM 2006 55(6):1219-26.

## Usage

``` r
scale_amp_molal(
  fit_result,
  ref_data,
  te,
  tr,
  water_t1,
  water_t2,
  metab_t1,
  metab_t2,
  ...
)
```

## Arguments

- fit_result:

  result object generated from fitting.

- ref_data:

  water reference MRS data object.

- te:

  the MRS TE in seconds.

- tr:

  the MRS TR in seconds.

- water_t1:

  assumed water T1 value.

- water_t2:

  assumed water T2 value.

- metab_t1:

  assumed metabolite T1 value.

- metab_t2:

  assumed metabolite T2 value.

- ...:

  additional arguments to get_td_amp function.

## Value

A `fit_result` object with a rescaled results table.
