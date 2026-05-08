# Make a basis-set object from a directory containing LCModel formatted RAW files.

Make a basis-set object from a directory containing LCModel formatted
RAW files.

## Usage

``` r
make_basis_from_raw(dir_path, ft, fs, ref)
```

## Arguments

- dir_path:

  path to the directory containing LCModel RAW files. One file per
  signal.

- ft:

  transmitter frequency in Hz.

- fs:

  sampling frequency in Hz.

- ref:

  reference value for ppm scale.

## Value

a basis-set object.
