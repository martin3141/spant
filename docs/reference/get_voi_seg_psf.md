# Return the white matter, gray matter and CSF composition of a volume.

Return the white matter, gray matter and CSF composition of a volume.

## Usage

``` r
get_voi_seg_psf(psf, mri_seg)
```

## Arguments

- psf:

  volume data as a nifti object.

- mri_seg:

  segmented brain volume as a nifti object.

## Value

a vector of partial volumes expressed as percentages.
