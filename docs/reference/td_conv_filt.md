# Time-domain convolution based filter.

Time-domain convolution based filter described by: Marion D, Ikura M,
Bax A. Improved solvent suppression in one-dimensional and
twodimensional NMR spectra by convolution of time-domain data. J Magn
Reson 1989;84:425-430.

## Usage

``` r
td_conv_filt(mrs_data, K = 25, ext = 1)
```

## Arguments

- mrs_data:

  MRS data to be filtered.

- K:

  window width in data points.

- ext:

  point separation for linear extrapolation.
