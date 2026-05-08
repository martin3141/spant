# Plot a volume as an overlay on a segmented brain volume.

Plot a volume as an overlay on a segmented brain volume.

## Usage

``` r
plot_voi_overlay_seg(mri_seg, voi, export_path = NULL, ...)
```

## Arguments

- mri_seg:

  segmented brain volume as a nifti object.

- voi:

  volume data as a nifti object.

- export_path:

  optional path to save the image in png format.

- ...:

  additional arguments to the ortho3 function.
