# Apply line-broadening (apodisation) to MRS data or basis object.

Apply line-broadening (apodisation) to MRS data or basis object.

## Usage

``` r
lb(x, lb, lg = 1)

# S3 method for class 'list'
lb(x, lb, lg = 1)

# S3 method for class 'mrs_data'
lb(x, lb, lg = 1)

# S3 method for class 'basis_set'
lb(x, lb, lg = 1)
```

## Arguments

- x:

  input mrs_data or basis_set object.

- lb:

  amount of line-broadening in Hz.

- lg:

  Lorentz-Gauss lineshape parameter (between 0 and 1).

## Value

line-broadened data.
