# Decimate an MRS signal by filtering in the time domain before downsampling.

Decimate an MRS signal by filtering in the time domain before
downsampling.

## Usage

``` r
decimate_mrs_td(mrs_data, q = 2, n = 4, ftype = "iir")
```

## Arguments

- mrs_data:

  MRS data object.

- q:

  integer factor to downsample by (default = 2).

- n:

  filter order used in the downsampling.

- ftype:

  filter type, "iir" or "fir".

## Value

decimated data.
