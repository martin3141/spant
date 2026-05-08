# Segment T1 weighted MRI data using FSL FAST and write to file. Runs bet as a preprocessing step by default.

This function requires a working installation of FSL and uses the fslr
package. You may need to specify the fsl install directory, eg:
'options(fsl.path = "/path/to/fsl")'

## Usage

``` r
segment_t1_fsl(mri_path, out_dir = NULL, deface = FALSE, bet_fit = 0.5)
```

## Arguments

- mri_path:

  path to the volumetric T1 data.

- out_dir:

  optional output directory. Defaults to the same directory as mri_path
  if not specified.

- deface:

  deface the input T1 data before analysis. Defaults to FALSE.

- bet_fit:

  fractional intensity threshold for bet brain extraction. Values should
  be between 0 and 1. Defaults to 0.5 with smaller values giving larger
  brain estimates.
