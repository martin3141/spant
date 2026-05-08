# Split paths into parts based on backslash, forwardslash and dot characters and return a data frame of these parts.

Split paths into parts based on backslash, forwardslash and dot
characters and return a data frame of these parts.

## Usage

``` r
paths2df(paths, col_num = NULL, extra_regex = NULL)
```

## Arguments

- paths:

  a vector of paths.

- col_num:

  when set, only the specified column number of the data frame will be
  returned.

- extra_regex:

  extra regular expression for splitting the paths. Example, "\_\|-"
  could be used to additionally split on underscores and hyphen
  characters.

## Value

data.frame of separated path elements.
