# Display an interactive orthographic projection plot of a nifti object.

Display an interactive orthographic projection plot of a nifti object.

## Usage

``` r
ortho3_inter(
  underlay,
  overlay = NULL,
  xyz = NULL,
  zlim = NULL,
  zlim_ol = NULL,
  alpha = 0.7,
  ...
)
```

## Arguments

- underlay:

  underlay image to be shown in grayscale.

- overlay:

  optional overlay image.

- xyz:

  x, y, z slice coordinates to display.

- zlim:

  underlay intensity limits.

- zlim_ol:

  overlay intensity limits.

- alpha:

  transparency of overlay.

- ...:

  other options to be passed to the ortho3 function.
