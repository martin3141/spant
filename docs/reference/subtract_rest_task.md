# Subtract mean rest spectrum from mean task spectrum after applying optimal linebroadening to the mean task spectrum. Usually used to correct for the BOLD lineshape narrowing effect in 1H fMRS data.

Subtract mean rest spectrum from mean task spectrum after applying
optimal linebroadening to the mean task spectrum. Usually used to
correct for the BOLD lineshape narrowing effect in 1H fMRS data.

## Usage

``` r
subtract_rest_task(mrs_data, task_vec, xlim = c(2.11, 1.91))
```

## Arguments

- mrs_data:

  dynamic MRS data.

- task_vec:

  a logical vector with the same length and dynamic scans. Elements set
  to TRUE and FALSE will be assigned to task and rest respectively.

- xlim:

  spectral region to match, eg c(2.11, 1.91) for tNAA region.

## Value

a list containing the matched mrs_data, difference spectra and
optimisation results.
