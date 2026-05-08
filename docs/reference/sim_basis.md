# Simulate a basis set object.

Simulate a basis set object.

## Usage

``` r
sim_basis(
  mol_list,
  pul_seq = seq_pulse_acquire,
  acq_paras = def_acq_paras(),
  xlim = NULL,
  auto_scale = FALSE,
  use_basis_cache = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- mol_list:

  list of `mol_parameter` objects. Alternatively, a character vector
  matching molecules may also be provided. Use the get_mol_names
  function for a full list of molecules.

- pul_seq:

  pulse sequence function to use.

- acq_paras:

  list of acquisition parameters or an mrs_data object. See
  [`def_acq_paras`](https://martin3141.github.io/spant/reference/def_acq_paras.md)

- xlim:

  ppm range limiting signals to be simulated.

- auto_scale:

  scale the basis based on the intensity of a singlet resonance. Needed
  for sequences with spatial simulation.

- use_basis_cache:

  create and use a cache of simulated basis sets stored in the
  "basis_cache" folder in the spant resources directory. Defaults to
  FALSE.

- verbose:

  output simulation progress and timings.

- ...:

  extra parameters to pass to the pulse sequence function.

## Value

basis object.
