# Trim a vector of filesystem paths.

Trim a vector of filesystem paths.

## Usage

``` r
trim_paths(paths, dir = 0, char = 0)
```

## Arguments

- paths:

  vectors of filesystem paths.

- dir:

  number of times to apply the base dirname function.

- char:

  number of characters to trim from the end of each path. Note this is
  performed after the dir based trimming.

## Value

a vector of trimmed paths.
