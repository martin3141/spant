# Set the masked voxels in a 2D MRSI dataset to given spectrum.

Set the masked voxels in a 2D MRSI dataset to given spectrum.

## Usage

``` r
set_mask_xy_mat(mrs_data, mask, mask_mrs_data)
```

## Arguments

- mrs_data:

  MRSI data object.

- mask:

  matrix of boolean values specifying the voxels to set, where a value
  of TRUE indicates the voxel should be set to mask_mrs_data.

- mask_mrs_data:

  the spectral data to be assigned to the masked voxels.

## Value

updated dataset.
