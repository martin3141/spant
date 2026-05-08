# Return a list of options for VARPRO based fitting.

Return a list of options for VARPRO based fitting.

## Usage

``` r
varpro_opts(
  nstart = 20,
  init_g_damping = 2,
  maxiters = 200,
  max_shift = 5,
  max_g_damping = 5,
  max_ind_damping = 5,
  anal_jac = TRUE,
  bl_smth_pts = 80
)
```

## Arguments

- nstart:

  position in the time-domain to start fitting, units of data points.

- init_g_damping:

  starting value for the global Gaussian line-broadening term - measured
  in Hz.

- maxiters:

  maximum number of levmar iterations to perform.

- max_shift:

  maximum shift allowed to each element in the basis set, measured in
  Hz.

- max_g_damping:

  maximum permitted global Gaussian line-broadening.

- max_ind_damping:

  maximum permitted Lorentzian line-broadening for each element in the
  basis set, measured in Hz.

- anal_jac:

  option to use the analytic or numerical Jacobian (logical).

- bl_smth_pts:

  number of data points to use in the baseline smoothing calculation.

## Value

list of options.

## Examples

``` r
varpro_opts(nstart = 10)
#> $nstart
#> [1] 10
#> 
#> $init_g_damping
#> [1] 2
#> 
#> $maxiters
#> [1] 200
#> 
#> $max_shift
#> [1] 5
#> 
#> $max_g_damping
#> [1] 5
#> 
#> $max_ind_damping
#> [1] 5
#> 
#> $anal_jac
#> [1] TRUE
#> 
#> $bl_smth_pts
#> [1] 80
#> 
```
