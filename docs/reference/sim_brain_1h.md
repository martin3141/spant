# Simulate MRS data with a similar appearance to normal brain (by default).

Simulate MRS data with a similar appearance to normal brain (by
default).

## Usage

``` r
sim_brain_1h(
  acq_paras = def_acq_paras(),
  type = "normal_v2",
  pul_seq = seq_slaser_ideal,
  xlim = c(0.5, 4.2),
  full_output = FALSE,
  amps = NULL,
  basis_lb = NULL,
  zero_lip_mm = FALSE,
  remove_lip_mm = FALSE,
  ...
)
```

## Arguments

- acq_paras:

  list of acquisition parameters or an mrs_data object. See
  [`def_acq_paras`](https://martin3141.github.io/spant/reference/def_acq_paras.md).

- type:

  type of spectrum, only "normal" is implemented currently.

- pul_seq:

  pulse sequence function to use.

- xlim:

  range of frequencies to simulate in ppm.

- full_output:

  when FALSE (default) only output the simulated MRS data. When TRUE
  output a list containing the MRS data, basis set object and
  corresponding amplitudes.

- amps:

  a vector of basis amplitudes may be specified to modify the output
  spectrum.

- basis_lb:

  apply additional Gaussian line-broadening to the basis (Hz).

- zero_lip_mm:

  zero the amplitudes of any lipid or macromolecular components based on
  their name starting with "MM" or "Lip".

- remove_lip_mm:

  remove any lipid or macromolecular basis components based on their
  name starting with "MM" or "Lip".

- ...:

  extra parameters to pass to the pulse sequence function.

## Value

see full_output option.
