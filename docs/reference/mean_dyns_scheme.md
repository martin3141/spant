# Average sets of dynamics according to a scheme vector.

Average sets of dynamics according to a scheme vector.

## Usage

``` r
mean_dyns_scheme(mrs_data, av_scheme)
```

## Arguments

- mrs_data:

  dynamic MRS data.

- av_scheme:

  vector containing consecutive integer values (starting at 1)
  representing sets of dynamics to average. This vector must have the
  same length as the number of dynamic scans in mrs_data. Any elements
  set to NA will not contribute to averaging.

## Value

Dynamic MRS data set containing the averaged dynamic sets.
