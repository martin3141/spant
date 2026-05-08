# Apodise MRSI data in the x-y direction with a k-space filter.

Apodise MRSI data in the x-y direction with a k-space filter.

## Usage

``` r
apodise_xy(mrs_data, func = "hamming", w = 2.5)
```

## Arguments

- mrs_data:

  MRSI data.

- func:

  must be "hamming", "hanning" or "gaussian".

- w:

  the reciprocal of the standard deviation for the Gaussian function.

## Value

apodised data.
