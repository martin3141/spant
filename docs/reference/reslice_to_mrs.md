# Reslice a nifti object to match the orientation of mrs data.

Reslice a nifti object to match the orientation of mrs data.

## Usage

``` r
reslice_to_mrs(mri, mrs, interp = 3L)
```

## Arguments

- mri:

  nifti object to be resliced.

- mrs:

  mrs_data object for the target orientation.

- interp:

  interpolation parameter, see nifyreg.linear definition.

## Value

resliced imaging data.
