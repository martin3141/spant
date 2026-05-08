# Perform a GLM analysis of dynamic MRS data in the spectral domain.

Perform a GLM analysis of dynamic MRS data in the spectral domain.

## Usage

``` r
glm_spec(mrs_data, regressor_df, full_output = FALSE)
```

## Arguments

- mrs_data:

  single-voxel dynamics MRS data.

- regressor_df:

  a data frame containing temporal regressors to be applied to each
  spectral datapoint.

- full_output:

  append mrs_data and regressor_df to the output list.

## Value

list of statistical results.
