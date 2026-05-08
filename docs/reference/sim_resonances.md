# Simulate a MRS data object containing a set of simulated resonances.

Simulate a MRS data object containing a set of simulated resonances.

## Usage

``` r
sim_resonances(
  freq = 0,
  amp = 1,
  lw = 0,
  lg = 0,
  phase = 0,
  freq_ppm = TRUE,
  acq_paras = def_acq_paras(),
  fp_scale = TRUE,
  back_extrap_pts = 0,
  sum_resonances = TRUE
)
```

## Arguments

- freq:

  resonance frequency.

- amp:

  resonance amplitude.

- lw:

  line width in Hz.

- lg:

  Lorentz-Gauss lineshape parameter (between 0 and 1).

- phase:

  phase in degrees.

- freq_ppm:

  frequencies are given in ppm units if set to TRUE, otherwise Hz are
  assumed.

- acq_paras:

  list of acquisition parameters. See
  [`def_acq_paras`](https://martin3141.github.io/spant/reference/def_acq_paras.md)

- fp_scale:

  multiply the first data point by 0.5.

- back_extrap_pts:

  number of data points to back extrapolate.

- sum_resonances:

  sum all resonances (default is TRUE), otherwise return a dynamic
  mrs_data object.

## Value

MRS data object.

## Examples

``` r
sim_data <- sim_resonances(freq = 2, lw = 5)
```
