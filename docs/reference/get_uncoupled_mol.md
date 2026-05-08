# Generate a `mol_parameters` object for a simple spin system with one resonance.

Generate a `mol_parameters` object for a simple spin system with one
resonance.

## Usage

``` r
get_uncoupled_mol(
  name,
  chem_shift,
  nucleus,
  scale_factor,
  lw,
  lg,
  full_name = NULL
)
```

## Arguments

- name:

  abbreviated name of the molecule.

- chem_shift:

  chemical shift of the resonance (PPM).

- nucleus:

  nucleus (1H, 31P...).

- scale_factor:

  multiplicative scaling factor. Note, this value can be made complex to
  adjust the phase of the resonance.

- lw:

  linewidth in Hz.

- lg:

  Lorentz-Gauss lineshape parameter (between 0 and 1).

- full_name:

  long name of the molecule (optional).

## Value

mol_parameters object.
