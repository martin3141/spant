# Generate a MRSI voxel PSF from an `mrs_data` object.

Generate a MRSI voxel PSF from an `mrs_data` object.

## Usage

``` r
get_mrsi_voxel_xy_psf(mrs_data, target_mri, x_pos, y_pos, z_pos)
```

## Arguments

- mrs_data:

  MRS data.

- target_mri:

  optional image data to match the intended volume space.

- x_pos:

  x voxel coordinate.

- y_pos:

  y voxel coordinate.

- z_pos:

  z voxel coordinate.

## Value

volume data as a nifti object.
