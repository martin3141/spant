# Display an orthographic projection plot of a nifti object.

Display an orthographic projection plot of a nifti object.

## Usage

``` r
ortho3(
  underlay,
  overlay = NULL,
  xyz = NULL,
  zlim = NULL,
  zlim_ol = NULL,
  alpha = 0.7,
  col_ol = viridisLite::viridis(64),
  orient_lab = TRUE,
  rescale = 1,
  crosshairs = TRUE,
  ch_lwd = 1,
  colourbar = TRUE,
  bg = "black",
  mar = c(0, 0, 0, 0),
  smallplot = c(0.63, 0.65, 0.07, 0.42),
  legend_axis_cex = 0.75
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

- col_ol:

  colour palette of overlay.

- orient_lab:

  display orientation labels (default TRUE).

- rescale:

  rescale factor for the underlay and overlay images.

- crosshairs:

  display the crosshairs (default TRUE).

- ch_lwd:

  crosshair linewidth.

- colourbar:

  display a colourbar for the overlay (default TRUE).

- bg:

  plot background colour.

- mar:

  plot margins.

- smallplot:

  smallplot option for positioning the colourbar.

- legend_axis_cex:

  font expansion factor for the legend axis text.
