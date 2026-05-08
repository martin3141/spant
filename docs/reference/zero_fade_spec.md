# Fade a spectrum to zero by frequency domain multiplication with a tanh function. Note this operation distorts data points at the end of the FID.

Fade a spectrum to zero by frequency domain multiplication with a tanh
function. Note this operation distorts data points at the end of the
FID.

## Usage

``` r
zero_fade_spec(mrs_data, start_ppm, end_ppm)
```

## Arguments

- mrs_data:

  data to be faded.

- start_ppm:

  start point of the fade in ppm units.

- end_ppm:

  end point of the fade in ppm units.

## Value

modified mrs_data object.
