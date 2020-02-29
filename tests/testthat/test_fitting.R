context("fitting")

test_that("ABfit no-optim works", {
  fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
                       package = "spant")
  
  mrs_data <- read_mrs(fname, format = "spar_sdat")
  basis    <- sim_basis_1h_brain_press(mrs_data)
  
  # don't run any optimisation steps
  opts <- abfit_opts(maxiters_pre = 0, maxiters = 0)
  fit_res  <- fit_mrs(mrs_data, basis, method = "abfit", opts = opts,
                      time = FALSE)
  
  expect_equal_to_reference(fit_res, "fit_res_abfit_no_optim.rds",
                            tolerance = 1e-4)
})

test_that("ABfit pre-fitting works", {
  fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
                       package = "spant")
  
  mrs_data <- read_mrs(fname, format = "spar_sdat")
  basis    <- sim_basis_1h_brain_press(mrs_data)
  
  # don't run the fine fitting step
  opts <- abfit_opts(maxiters = 0)
  fit_res  <- fit_mrs(mrs_data, basis, method = "abfit", opts = opts,
                      time = FALSE)
  
  expect_equal_to_reference(fit_res, "fit_res_abfit_pre.rds", tolerance = 1e-4)
})

test_that("ABfit full works", {
  fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
                       package = "spant")
  
  mrs_data <- read_mrs(fname, format = "spar_sdat")
  basis    <- sim_basis_1h_brain_press(mrs_data)
  
  fit_res  <- fit_mrs(mrs_data, basis, method = "abfit", time = FALSE)
  expect_equal_to_reference(fit_res, "fit_res_abfit.rds", tolerance = 1e-4)
})