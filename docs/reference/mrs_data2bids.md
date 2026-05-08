# Create a BIDS file structure from a vector of MRS data paths or list of mrs_data objects.

Create a BIDS file structure from a vector of MRS data paths or list of
mrs_data objects.

## Usage

``` r
mrs_data2bids(
  mrs_data,
  output_dir,
  suffix = NULL,
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
  skip_existing = TRUE
)
```

## Arguments

- mrs_data:

  vector of MRS data paths or list of mrs_data objects.

- output_dir:

  the base directory to create the BIDS structure.

- suffix:

  optional vector of file suffixes. Default behaviour is to
  automatically determine these from the input data, however it is
  recommended that they are specified to allow more efficient skipping
  of existing data.

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
