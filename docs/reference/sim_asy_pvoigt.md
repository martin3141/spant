# Generate an asymmetric pseudo-Voigt resonance in the frequency domain.

Method is described in detail by Stancik AL and Brauns EB: "A simple
asymmetric lineshape for fitting infrared absorption spectra." Vib
Spectrosc. 2008; 47: 66-69.

## Usage

``` r
sim_asy_pvoigt(
  freq = 0,
  fwhm = 0,
  lg = 0,
  asy = 0,
  acq_paras = def_acq_paras(),
  gen_im_pts = FALSE
)
```

## Arguments

- freq:

  resonance frequency in ppm.

- fwhm:

  resonance FWHM in Hz.

- lg:

  Lorentz-Gauss lineshape parameter (between 0 and 1).

- asy:

  asymmetry parameter.

- acq_paras:

  list of acquisition parameters. See
  [`def_acq_paras`](https://martin3141.github.io/spant/reference/def_acq_paras.md)

- gen_im_pts:

  option to generate imaginary data points, defaults to FALSE.
