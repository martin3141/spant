# Convert SPM style segmentation files to a single categorical image where the numerical values map as: 0) Other, 1) CSF, 2) GM and 3) WM.

Convert SPM style segmentation files to a single categorical image where
the numerical values map as: 0) Other, 1) CSF, 2) GM and 3) WM.

## Usage

``` r
spm_pve2categorical(fname)
```

## Arguments

- fname:

  any of the segmentation files (eg c1_MY_T1.nii).

## Value

nifti object.
