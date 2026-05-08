# Mask an MRS dataset in the dynamic dimension.

Mask an MRS dataset in the dynamic dimension.

## Usage

``` r
mask_dyns(mrs_data, mask)
```

## Arguments

- mrs_data:

  MRS data object.

- mask:

  vector of boolean values specifying the dynamics to mask, where a
  value of TRUE indicates the spectrum should be removed.

## Value

masked dataset.
