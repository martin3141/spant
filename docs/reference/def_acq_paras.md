# Return (and optionally modify using the input arguments) a list of the default acquisition parameters.

Return (and optionally modify using the input arguments) a list of the
default acquisition parameters.

## Usage

``` r
def_acq_paras(
  ft = getOption("spant.def_ft"),
  fs = getOption("spant.def_fs"),
  N = getOption("spant.def_N"),
  ref = getOption("spant.def_ref"),
  nuc = getOption("spant.def_nuc")
)
```

## Arguments

- ft:

  specify the transmitter frequency in Hz.

- fs:

  specify the sampling frequency in Hz.

- N:

  specify the number of data points in the spectral dimension.

- ref:

  specify the reference value for ppm scale.

- nuc:

  specify the resonant nucleus.

## Value

A list containing the following elements:

- ft transmitter frequency in Hz.

- fs sampling frequency in Hz.

- N number of data points in the spectral dimension.

- ref reference value for ppm scale.

- nuc resonant nucleus.
