# Read MRS data from the filesystem.

Read MRS data from the filesystem.

## Usage

``` r
read_mrs(
  path,
  format = NULL,
  ft = NULL,
  fs = NULL,
  ref = NULL,
  n_ref_scans = NULL,
  full_fid = FALSE,
  omit_svs_ref_scans = TRUE,
  verbose = FALSE,
  extra = NULL,
  fid_filt_dist = NULL
)
```

## Arguments

- path:

  file name or directory containing the MRS data.

- format:

  string describing the data format. Must be one of the following :
  "spar_sdat", "rda", "dicom", "twix", "pfile", "list_data", "paravis",
  "dpt", "lcm_raw", "rds", "nifti", "varian", "jmrui_txt". If not
  specified, the format will be guessed from the filename extension, or
  will be assumed to be a Siemens ima dynamic data if the path is a
  directory.

- ft:

  transmitter frequency in Hz (required for list_data format).

- fs:

  sampling frequency in Hz (required for list_data format).

- ref:

  reference value for ppm scale (required for list_data format).

- n_ref_scans:

  override the number of water reference scans detected in the file
  header (GE p-file only).

- full_fid:

  export all data points, including those before the start of the FID
  (default = FALSE), TWIX format only.

- omit_svs_ref_scans:

  remove any reference scans sometimes saved in SVS twix data (default =
  TRUE).

- verbose:

  print data file information (default = FALSE).

- extra:

  an optional data frame to provide additional variables for use in
  subsequent analysis steps, eg id or grouping variables.

- fid_filt_dist:

  indicate if the data has a distorted FID due to a brick-wall filter
  being used to downsample the data. Default is to auto detect this from
  the data, but TRUE or FALSE options can be given to override
  detection.

## Value

MRS data object.

## Examples

``` r
fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")
mrs_data <- read_mrs(fname)
print(mrs_data)
#> MRS Data Parameters
#> ----------------------------------
#> Trans. freq (MHz)       : 127.7861
#> FID data points         : 1024
#> X,Y,Z dimensions        : 1x1x1
#> Dynamics                : 1
#> Coils                   : 1
#> Voxel resolution (mm)   : 20x20x20
#> Sampling frequency (Hz) : 2000
#> Repetition time (s)     : 2 
#> Reference freq. (ppm)   : 4.65
#> Nucleus                 : 1H
#> Spectral domain         : FALSE
#> Number of transients    : 128 
#> Echo time (s)           : 0.03 
#> Manufacturer            : Philips 
```
