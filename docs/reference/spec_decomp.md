# Decompose an mrs_data object into white and gray matter spectra.

An implementation of the method published by Goryawala et al MRM 79(6)
2886-2895 (2018). "Spectral decomposition for resolving partial volume
effects in MRSI".

## Usage

``` r
spec_decomp(mrs_data, wm, gm, norm_fractions = TRUE)
```

## Arguments

- mrs_data:

  data to be decomposed into white and gray matter spectra.

- wm:

  vector of white matter contributions to each voxel.

- gm:

  vector of gray matter contributions to each voxel.

- norm_fractions:

  option to normalise the wm, gm vectors for each voxel.

## Value

a list of two mrs_data objects corresponding to the two tissue types.
