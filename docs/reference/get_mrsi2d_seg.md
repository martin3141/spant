# Calculate the partial volume estimates for each voxel in a 2D MRSI dataset.

Localisation is assumed to be perfect in the z direction and determined
by the ker input in the x-y direction.

## Usage

``` r
get_mrsi2d_seg(mrs_data, mri_seg, ker)
```

## Arguments

- mrs_data:

  2D MRSI data with multiple voxels in the x-y dimension.

- mri_seg:

  MRI data with values corresponding to the segmentation class. Must be
  1mm isotropic resolution.

- ker:

  MRSI PSF kernel in the x-y direction compatible with the mmand
  package, eg: mmand::shapeKernel(c(10, 10), type = "box").

## Value

a data frame of partial volume estimates and individual segmentation
maps.
