# Apply water reference scaling to a fitting results object to yield metabolite quantities in units of "mmol per Kg wet weight".

See the LCModel manual (section 10.2) on water-scaling for details on
the assumptions and relevant references. Use this type of concentration
scaling to compare fit results with LCModel and TARQUIN defaults.
Otherwise scale_amp_molal_pvc is the preferred method. Note, the LCModel
manual (section 1.3) states:

## Usage

``` r
scale_amp_legacy(fit_result, ref_data, w_att = 0.7, w_conc = 35880, ...)
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

## Details

"Concentrations should be labelled 'mmol per Kg wet weight'. We use the
shorter (incorrect) abbreviation mM. The actual mM is the mmol per Kg
wet weight multiplied by the specific gravity of the tissue, typically
1.04 in brain."
