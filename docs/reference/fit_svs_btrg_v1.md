# A simple wrapper of fit_svs with defaults adjusted to improve the analysis of lower quality brain tumour spectra.

A simple wrapper of fit_svs with defaults adjusted to improve the
analysis of lower quality brain tumour spectra.

## Usage

``` r
fit_svs_btrg_v1(
  append_basis = c("peth", "cit", "gly"),
  pre_align_ref_freq = c(4.65),
  hsvd_width = 50,
  fit_opts = NULL,
  ...
)
```

## Arguments

- append_basis:

  names of extra signals to add to the default basis. Eg append_basis =
  c("peth", "cit"). Use get_mol_names() function to print all available
  signals. Cannot be used with precompiled basis sets.

- pre_align_ref_freq:

  reference frequency in ppm units. More than one frequency may be
  specified. Defaults to : c(4.65).

- hsvd_width:

  set the width of the HSVD filter in Hz. Note the applied width is
  between -width and +width Hz, with 0 Hz being defined at the centre of
  the spectral width. Default is set to 50 Hz.

- fit_opts:

  options to pass to the fitting method.

- ...:

  other arguments to pass to fit_svs.
