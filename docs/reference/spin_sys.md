# Create a spin system object for pulse sequence simulation.

Create a spin system object for pulse sequence simulation.

## Usage

``` r
spin_sys(spin_params, ft, ref, precomp_jc_H = NULL, precomp_Iz = NULL)
```

## Arguments

- spin_params:

  an object describing the spin system properties.

- ft:

  transmitter frequency in Hz.

- ref:

  reference value for ppm scale.

- precomp_jc_H:

  use a precomputed J-coupling H matrix to save time.

- precomp_Iz:

  use precomputed Iz matrices to save time.

## Value

spin system object.
