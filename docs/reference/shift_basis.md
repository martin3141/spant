# Apply frequency shifts to basis set signals.

Apply frequency shifts to basis set signals.

## Usage

``` r
shift_basis(basis, shifts)
```

## Arguments

- basis:

  the basis to apply the shift to.

- shifts:

  a vector of frequency shifts to apply in ppm units. Must be the same
  length as there are basis elements, or one value to be applied to all
  elements.

## Value

modified basis set object.
