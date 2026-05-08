# Papoulis-Gerchberg (PG) algorithm method for k-space extrapolation.

PG method as described in: Haupt CI, Schuff N, Weiner MW, Maudsley AA.
Removal of lipid artifacts in 1H spectroscopic imaging by data
extrapolation. Magn Reson Med. 1996 May;35(5):678-87. Extrapolation is
performed to expand k-space coverage by a factor of 2, with the aim to
reduce Gibbs ringing.

## Usage

``` r
pg_extrap_xy(
  mrs_data,
  img_mask = NULL,
  kspace_mask = NULL,
  intensity_thresh = 0.15,
  iters = 50
)
```

## Arguments

- mrs_data:

  MRS data object.

- img_mask:

  a boolean matrix of voxels with strong signals to be extrapolated.
  Must be twice the dimensions of the input data.

- kspace_mask:

  a boolean matrix of kspace points that have been sampled. Typically a
  circle for MRSI, but defaults to the full rectangular area of k-space
  covered by the input data. Must match the x-y dimensions of the input
  data.

- intensity_thresh:

  used to define img_mask based on the strength of the signal in each
  voxel. Defaults to intensities greater than 15% of the maximum.
  Ignored if img_mask is specified as argument.

- iters:

  number of iterations to perform.

## Value

extrapolated `mrs_data` object.
