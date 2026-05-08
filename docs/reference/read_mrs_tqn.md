# Read MRS data using the TARQUIN software package.

Read MRS data using the TARQUIN software package.

## Usage

``` r
read_mrs_tqn(fname, fname_ref = NA, format, id = NA, group = NA)
```

## Arguments

- fname:

  the filename containing the MRS data.

- fname_ref:

  a second filename containing reference MRS data.

- format:

  format of the MRS data. Can be one of the following: siemens, philips,
  ge, dcm, dpt, rda, lcm, varian, bruker, jmrui_txt.

- id:

  optional ID string.

- group:

  optional group string.

## Value

MRS data object.

## Examples

``` r
fname <- system.file("extdata","philips_spar_sdat_WS.SDAT",package="spant")
if (FALSE) { # \dontrun{
mrs_data <- read_mrs_tqn(fname, format="philips")
} # }
```
