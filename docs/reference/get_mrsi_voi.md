# Generate a MRSI VOI from an `mrs_data` object.

Generate a MRSI VOI from an `mrs_data` object.

## Usage

``` r
get_mrsi_voi(mrs_data, target_mri = NULL, map = NULL, ker = mmand::boxKernel())
```

## Arguments

- mrs_data:

  MRS data.

- target_mri:

  optional image data to match the intended volume space.

- map:

  optional voi intensity map.

- ker:

  kernel to rescale the map data to the target_mri. Default value is
  mmand::boxKernel(), use mmand::mnKernel() for a smoothed map.

## Value

volume data as a nifti object.
