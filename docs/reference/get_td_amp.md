# Return an array of amplitudes derived from fitting the initial points in the time domain and extrapolating back to t=0.

Return an array of amplitudes derived from fitting the initial points in
the time domain and extrapolating back to t=0.

## Usage

``` r
get_td_amp(mrs_data, nstart = 10, nend = 50, method = "poly")
```

## Arguments

- mrs_data:

  MRS data.

- nstart:

  first data point to fit.

- nend:

  last data point to fit.

- method:

  method for measuring the amplitude, one of "poly", spline" or exp".

## Value

array of amplitudes.
