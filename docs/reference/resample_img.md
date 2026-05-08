# Resample an image to match a target image space.

Resample an image to match a target image space.

## Usage

``` r
resample_img(source, target, interp = 3L)
```

## Arguments

- source:

  image data as a nifti object.

- target:

  image data as a nifti object.

- interp:

  interpolation parameter, see nifyreg.linear definition.

## Value

resampled image data as a nifti object.
