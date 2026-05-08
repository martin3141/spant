# Scale a basis-set to be consistent with spant assumptions for water scaling.

For correct water scaling, spant assumes the time-domain amplitude (t =
0) for a single proton is 0.5. Internally simulated basis-sets will be
correctly scaled, however imported basis-sets should be assumed to be
un-scaled and this function should be used. Note that the singlet
specified should only contain one resonance, and that any additional
signals (eg TSP or residual water) will result in incorrect scaling.
Therefore, only simulated basis sets are appropriate for use with this
function.

## Usage

``` r
scale_basis_from_singlet(basis, name, protons)
```

## Arguments

- basis:

  basis set to be scaled.

- name:

  the name of the singlet to be used as a scaling reference.

- protons:

  the number of MRS visible protons contributing to the singlet
  resonance.

## Value

a scaled basis.
