# Create a BIDS file structure from a vector of data paths or list of mri/mrs data objects.

Create a BIDS file structure from a vector of data paths or list of
mri/mrs data objects.

## Usage

``` r
mr_data2bids(
  mr_data,
  suffix,
  output_dir,
  sub = NULL,
  ses = NULL,
  task = NULL,
  acq = NULL,
  nuc = NULL,
  voi = NULL,
  rec = NULL,
  run = NULL,
  echo = NULL,
  inv = NULL,
  skip_existing = TRUE,
  mri_format = "nifti",
  deface_mri = FALSE
)
```

## Arguments

- mr_data:

  vector of data paths or list of mri/mrs objects.

- suffix:

  vector of file suffixes, eg : c("svs", "mrsi", "T1w).

- output_dir:

  the base directory to create the BIDS structure.

- sub:

  optional vector of subject labels. If not specified, these will be
  automatically generated as a series of increasing zero-padded integer
  values corresponding to the mrs_data input indices.

- ses:

  optional vector of session labels.

- task:

  optional vector of task labels.

- acq:

  optional vector of acquisition labels.

- nuc:

  optional vector of nucleus labels.

- voi:

  optional vector of volume of interest labels.

- rec:

  optional vector of reconstruction labels.

- run:

  optional vector of run indices.

- echo:

  optional vector of echo time indices.

- inv:

  optional vector of inversion indices.

- skip_existing:

  skip any data files that have already been converted. Defaults to
  TRUE, set to FALSE to force an overwrite of any existing data files.

- mri_format:

  defaults to "nifti", can also be "dicom" provided the divest package
  is installed.

- deface_mri:

  option to apply fsl_deface to the mri as a preprocessing step.
  Defaults to FALSE, requires ANTs and faceoff to be installed when
  TRUE.
