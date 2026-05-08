# Set the default acquisition parameters.

Set the default acquisition parameters.

## Usage

``` r
set_def_acq_paras(
  ft = getOption("spant.def_ft"),
  fs = getOption("spant.def_fs"),
  N = getOption("spant.def_N"),
  ref = getOption("spant.def_ref"),
  nuc = getOption("spant.nuc")
)
```

## Arguments

- ft:

  transmitter frequency in Hz.

- fs:

  sampling frequency in Hz.

- N:

  number of data points in the spectral dimension.

- ref:

  reference value for ppm scale.

- nuc:

  resonant nucleus.
