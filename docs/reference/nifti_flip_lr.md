# Flip the x data dimension order of a nifti image. This corresponds to flipping MRI data in the left-right direction, assuming the data in save in neurological format (can check with fslorient program).

Flip the x data dimension order of a nifti image. This corresponds to
flipping MRI data in the left-right direction, assuming the data in save
in neurological format (can check with fslorient program).

## Usage

``` r
nifti_flip_lr(x)
```

## Arguments

- x:

  nifti object to be processed.

## Value

nifti object with reversed x data direction.
