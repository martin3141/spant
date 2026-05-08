# Combine fitting results for group analysis.

Combine fitting results for group analysis.

## Usage

``` r
fit_svs_group_results(
  search_path = NULL,
  paths = NULL,
  output_dir = "fit_svs_group_results",
  verbose = TRUE
)
```

## Arguments

- search_path:

  path to start recursive search for fitting results. Cannot be used
  together with the paths argument.

- paths:

  a set of paths to spant output files, usually named :
  "spant_fit_svs_data.rds". Cannot be used together with the search_path
  argument.

- output_dir:

  directory path to store group results.

- verbose:

  verbose, defaults to TRUE.
