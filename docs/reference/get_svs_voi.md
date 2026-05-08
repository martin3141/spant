# Generate a SVS acquisition volume from an `mrs_data` object.

Generate a SVS acquisition volume from an `mrs_data` object.

## Usage

``` r
get_svs_voi(mrs_data, target_mri)
```

## Arguments

- mrs_data:

  MRS data.

- target_mri:

  optional image data to match the intended volume space.

## Value

volume data as a nifti object.
