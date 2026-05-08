# Read the text format header found in Siemens IMA and TWIX data files.

Read the text format header found in Siemens IMA and TWIX data files.

## Usage

``` r
read_siemens_txt_hdr(input, version = "vd", verbose = FALSE, offset = 0)
```

## Arguments

- input:

  file name to read or raw data.

- version:

  software version, can be "vb" or "vd".

- verbose:

  print information to the console.

- offset:

  offset to begin searching for the text header.

## Value

a list of parameter values
