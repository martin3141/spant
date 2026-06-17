# Remove elements from a basis set object.

Remove elements from a basis set object.

## Usage

``` r
rm_basis_elements(basis, rm_str)
```

## Arguments

- basis:

  input basis.

- rm_str:

  a grep expression to remove basis elements. Use "\|" for multiple
  matches, eg to match alanine and lactate only : "^Ala\$\|^Lac\$".

## Value

basis with elements removed.
