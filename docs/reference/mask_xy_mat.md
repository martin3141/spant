# Mask a 2D MRSI dataset in the x-y dimension.

Mask a 2D MRSI dataset in the x-y dimension.

## Usage

``` r
mask_xy_mat(mrs_data, mask, value = NA)
```

## Arguments

- mrs_data:

  MRS data object.

- mask:

  matrix of boolean values specifying the voxels to mask, where a value
  of TRUE indicates the voxel should be removed.

- value:

  the value to set masked data to (usually NA or 0).

## Value

masked dataset.
