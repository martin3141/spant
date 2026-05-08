# Generate a double gamma model of the HRF as used in SPM.

Generate a double gamma model of the HRF as used in SPM.

## Usage

``` r
get_hrf(end_t = 30, res_t = 0.01)
```

## Arguments

- end_t:

  last time point to generate in seconds.

- res_t:

  temporal resolution in seconds, defaults to 10ms.

## Value

a data.frame of time and HRF vectors.
