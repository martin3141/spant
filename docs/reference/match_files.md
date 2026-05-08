# Match files based on a vector of input paths and a glob pattern. The glob pattern is appended to each path and should match one file only.

Match files based on a vector of input paths and a glob pattern. The
glob pattern is appended to each path and should match one file only.

## Usage

``` r
match_files(paths, glob)
```

## Arguments

- paths:

  vectors of filesystem paths.

- glob:

  pattern to append to each path before passing to Sys.glob.

## Value

matched files. NA values correspond to either no match or multiple
matches.
