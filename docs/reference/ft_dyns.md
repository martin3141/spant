# Apply the Fourier transform over the dynamic dimension.

Apply the Fourier transform over the dynamic dimension.

## Usage

``` r
ft_dyns(mrs_data, ft_shift = FALSE, ret_mod = FALSE, fd = TRUE)
```

## Arguments

- mrs_data:

  MRS data where the dynamic dimension is in the time-domain.

- ft_shift:

  apply FT shift to the output, default is FALSE.

- ret_mod:

  return the modulus out the transform, default is FALSE.

- fd:

  transform the chemical shift axis to the frequency domain first,
  default is TRUE.

## Value

transformed MRS data.
