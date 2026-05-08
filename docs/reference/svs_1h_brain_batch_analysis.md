# Batch interface to the standard SVS 1H brain analysis pipeline.

Batch interface to the standard SVS 1H brain analysis pipeline.

## Usage

``` r
svs_1h_brain_batch_analysis(
  metab_list,
  w_ref_list = NULL,
  mri_seg_list = NULL,
  mri_list = NULL,
  output_dir_list = NULL,
  extra = NULL,
  ...
)
```

## Arguments

- metab_list:

  list of file paths or mrs_data objects containing MRS metabolite data.

- w_ref_list:

  list of file paths or mrs_data objects containing MRS water reference
  data.

- mri_seg_list:

  list of file paths or nifti objects containing segmented MRI data.

- mri_list:

  list of file paths or nifti objects containing anatomical MRI data.

- output_dir_list:

  list of directory paths to output fitting results.

- extra:

  a data frame with the same number of rows as metab_list, containing
  additional information to be attached to the fit results table.

- ...:

  additional options to be passed to the svs_1h_brain_analysis function.

## Value

a list of fit_result objects.
