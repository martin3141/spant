# Return the ANTs installation directory, or throw an error if not found.

Will check and return the "spant.ants_dir" option set by set_ants_dir.
If not set, will search the spant_resources directory for ANTs and
return the most recent version if multiple are found.

## Usage

``` r
get_ants_dir()
```

## Value

ANTs installation directory.
