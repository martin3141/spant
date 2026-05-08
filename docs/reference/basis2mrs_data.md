# Convert a basis object to an mrs_data object - where basis signals are spread across the dynamic dimension.

Convert a basis object to an mrs_data object - where basis signals are
spread across the dynamic dimension.

## Usage

``` r
basis2mrs_data(
  basis,
  sum_elements = FALSE,
  amps = NULL,
  shifts = NULL,
  lbs = NULL
)
```

## Arguments

- basis:

  basis set object.

- sum_elements:

  return the sum of basis elements (logical)

- amps:

  a vector of scaling factors to apply to each basis element.

- shifts:

  a vector of frequency shifts (in ppm) to apply to each basis element.

- lbs:

  a vector of Lorentzian line broadening terms (in Hz) to apply to each
  basis element.

## Value

an mrs_data object with basis signals spread across the dynamic
dimension or summed.
