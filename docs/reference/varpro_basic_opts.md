# Return a list of options for a basic VARPRO analysis.

Return a list of options for a basic VARPRO analysis.

## Usage

``` r
varpro_basic_opts(method = "fd_re", nnls = TRUE, ppm_left = 4, ppm_right = 0.2)
```

## Arguments

- method:

  one of "td", "fd", "fd_re".

- nnls:

  restrict basis amplitudes to non-negative values.

- ppm_left:

  downfield frequency limit for the fitting range (ppm).

- ppm_right:

  upfield frequency limit for the fitting range (ppm).

## Value

full list of options.
