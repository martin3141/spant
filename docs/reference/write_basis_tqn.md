# Generate a basis file using TARQUIN.

Generate a basis file using TARQUIN.

## Usage

``` r
write_basis_tqn(basis_file, metab_data, opts = NULL)
```

## Arguments

- basis_file:

  filename of the basis file to be generated.

- metab_data:

  MRS data object to match the generated basis parameters.

- opts:

  list of options to pass to TARQUIN.

## Examples

``` r
if (FALSE) { # \dontrun{
write_basis_tqn('test.basis',mrs_data,c("--echo","0.04"))
} # }
```
