# Deface a T1 weighted head scan using the FaceOff method described in : https://github.com/srikash/FaceOff. Requires ANTs to be installed.

Deface a T1 weighted head scan using the FaceOff method described in :
https://github.com/srikash/FaceOff. Requires ANTs to be installed.

## Usage

``` r
faceoff(mri_path, out_dir = NULL)
```

## Arguments

- mri_path:

  path to the volumetric T1 data.

- out_dir:

  optional output directory. Defaults to the same directory as mri_path
  if not specified.
