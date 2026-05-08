# Apply a function over specified array axes.

Apply a function over specified array axes.

## Usage

``` r
apply_axes(x, axes, fun, ...)
```

## Arguments

- x:

  an array.

- axes:

  a vector of axes to apply fun over.

- fun:

  function to be applied.

- ...:

  optional arguments to fun.

## Value

array.

## Examples

``` r
z <- array(1:1000, dim = c(10, 10, 10))
a <- apply_axes(z, 3, fft)
#> [1] 2 3 1
a[1,1,] == fft(z[1,1,])
#>  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
a <- apply_axes(z, 3, sum)
a[1,1,] == sum(z[1,1,])
#> [1] TRUE
```
