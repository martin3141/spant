# Plot a volume as an image overlay.

Plot a volume as an image overlay.

## Usage

``` r
plot_voi_overlay(mri, voi, export_path = NULL, zlim = NULL, ...)
```

## Arguments

- mri:

  image data as a nifti object or path to data file.

- voi:

  volume data as a nifti object or path to data file.

- export_path:

  optional path to save the image in png format.

- zlim:

  underlay intensity limits.

- ...:

  additional arguments to the ortho3 function.
