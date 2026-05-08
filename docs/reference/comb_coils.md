# Combine coil data based on the first data point of a reference signal.

By default, elements are phased and scaled prior to summation. Where a
reference signal is not given, the mean dynamic signal will be used
instead.

## Usage

``` r
comb_coils(
  metab,
  ref = NULL,
  noise = NULL,
  scale = TRUE,
  scale_method = "sig_noise_sq",
  sum_coils = TRUE,
  noise_region = c(-0.5, -2.5),
  average_ref_dyns = TRUE,
  ref_pt_index = 1,
  ret_metab_only = FALSE
)
```

## Arguments

- metab:

  MRS data containing metabolite data.

- ref:

  MRS data containing reference data (optional).

- noise:

  MRS data from a noise scan (optional).

- scale:

  option to rescale coil elements based on the first data point
  (logical).

- scale_method:

  one of "sig_noise_sq", "sig_noise" or "sig".

- sum_coils:

  sum the coil elements as a final step (logical).

- noise_region:

  the spectral region (in ppm) to estimate the noise.

- average_ref_dyns:

  take the mean of the reference scans in the dynamic dimension before
  use.

- ref_pt_index:

  time-domain point to use for estimating phase and scaling values.

- ret_metab_only:

  return the metabolite data only, even if reference data has been
  specified.

## Value

MRS data.
