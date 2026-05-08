# Perform a fit based analysis of MRS data.

Note that TARQUIN and LCModel require these packages to be installed,
and the functions set_tqn_cmd and set_lcm_cmd (respectively) need to be
used to specify the location of these software packages.

## Usage

``` r
fit_mrs(
  metab,
  basis = NULL,
  method = "ABFIT",
  w_ref = NULL,
  opts = NULL,
  parallel = FALSE,
  cl = NULL,
  time = TRUE,
  progress = "text",
  extra = NULL
)
```

## Arguments

- metab:

  metabolite data.

- basis:

  basis class object or character vector to basis file in LCModel .basis
  format.

- method:

  'ABFIT' (default), 'VARPRO', 'VARPRO_3P', 'TARQUIN' or 'LCMODEL'.

- w_ref:

  water reference data for concentration scaling (optional).

- opts:

  options to pass to the analysis method.

- parallel:

  perform analyses in parallel (TRUE or FALSE).

- cl:

  a parallel socket cluster required to run analyses in parallel. Eg, cl
  \<- parallel::makeCluster(4).

- time:

  measure the time taken for the analysis to complete (TRUE or FALSE).

- progress:

  option is passed to plyr::alply function to display a progress bar
  during fitting. Default value is "text", set to "none" to disable.

- extra:

  an optional data frame to provide additional variables for use in
  subsequent analysis steps, eg id or grouping variables.

## Value

MRS analysis object.

## Details

Fitting approaches described in the following references: ABfit Wilson,
M. Adaptive baseline fitting for 1H MR spectroscopy analysis. Magn Reson
Med 2012;85:13-29.

VARPRO van der Veen JW, de Beer R, Luyten PR, van Ormondt D. Accurate
quantification of in vivo 31P NMR signals using the variable projection
method and prior knowledge. Magn Reson Med 1988;6:92-98.

TARQUIN Wilson, M., Reynolds, G., Kauppinen, R. A., Arvanitis, T. N. &
Peet, A. C. A constrained least-squares approach to the automated
quantitation of in vivo 1H magnetic resonance spectroscopy data. Magn
Reson Med 2011;65:1-12.

LCModel Provencher SW. Estimation of metabolite concentrations from
localized in vivo proton NMR spectra. Magn Reson Med 1993;30:672-679.

## Examples

``` r
fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package =
"spant")
svs <- read_mrs(fname)
if (FALSE) { # \dontrun{
basis <- sim_basis_1h_brain_press(svs)
fit_result <- fit_mrs(svs, basis)
} # }
```
