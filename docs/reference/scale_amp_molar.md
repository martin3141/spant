# Apply water reference scaling to a fitting results object to yield metabolite quantities in millimolar (mM) units (mol / Litre of tissue). This function is depreciated, please use scale_amp_legacy instead.

See the LCModel manual (section 10.2) on water-scaling for details on
the assumptions and relevant references. Use this type of concentration
scaling to compare fit results with LCModel and TARQUIN defaults.
Otherwise scale_amp_molal_pvc is generally the preferred method.

## Usage

``` r
scale_amp_molar(fit_result, ref_data, w_att = 0.7, w_conc = 35880, ...)
```

## Arguments

- fit_result:

  a result object generated from fitting.

- ref_data:

  water reference MRS data object.

- w_att:

  water attenuation factor (default = 0.7). Assumes water T2 of 80ms and
  a TE = 30 ms. exp(-30ms / 80ms) ~ 0.7.

- w_conc:

  assumed water concentration (default = 35880). Default value
  corresponds to typical white matter. Set to 43300 for gray matter, and
  55556 for phantom measurements.

- ...:

  additional arguments to get_td_amp function.

## Value

a `fit_result` object with a rescaled results table.
