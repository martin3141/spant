# Back extrapolate time-domain data points using an autoregressive model.

Back extrapolate time-domain data points using an autoregressive model.

## Usage

``` r
back_extrap_ar(
  mrs_data,
  extrap_pts,
  pred_pts = NULL,
  method = "burg",
  rem_add = TRUE,
  ...
)
```

## Arguments

- mrs_data:

  mrs_data object.

- extrap_pts:

  number of points to extrapolate.

- pred_pts:

  number of points to base the extrapolation on.

- method:

  character string specifying the method to fit the model. Must be one
  of the strings in the default argument (the first few characters are
  sufficient). Defaults to "burg".

- rem_add:

  remove additional points from the end of the FID to maintain the
  original length of the dataset. Default to TRUE.

- ...:

  additional arguments to specific methods, see ?ar.

## Value

back extrapolated data.
