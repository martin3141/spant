context("fitting")

test_that("Test ABfit no-optim", {
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

test_that("Test ABfit coarse-fitting", {
  fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
                       package = "spant")
  
  mrs_data <- read_mrs(fname, format = "spar_sdat")
  basis    <- sim_basis_1h_brain_press(mrs_data)
  
  # don't run the fine fitting step
  opts <- abfit_opts(maxiters = 0)
  fit_res  <- fit_mrs(mrs_data, basis, method = "abfit", opts = opts,
                      time = FALSE)
  
  # this test is most sensitive to different platforms. Rounding errors can
  # cause different optimisation paths to be taken leading to differences in 
  # the number of iterations and small changes in the final solution. For these
  # reasons we only consider the minimum value found, which should be more
  # stable.
  expected <- 0.0001037292
  expect_equal(fit_res$res_tab$res.deviance, expected, tolerance = 1e-5,
               scale = expected)
})

test_that("Test ABfit fine-fitting", {
  fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
                       package = "spant")
  
  mrs_data <- read_mrs(fname, format = "spar_sdat")
  basis    <- sim_basis_1h_brain_press(mrs_data)
  
  # don't run the pre fitting step
  opts <- abfit_opts(maxiters_pre = 0)
  fit_res  <- fit_mrs(mrs_data, basis, method = "abfit", opts = opts,
                      time = FALSE)
  
  expect_equal_to_reference(fit_res, "fit_res_abfit_fine.rds", tolerance = 1e-4)
})

test_that("Test ABfit full", {
  fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
                       package = "spant")
  
  mrs_data <- read_mrs(fname, format = "spar_sdat")
  basis    <- sim_basis_1h_brain_press(mrs_data)
  
  fit_res  <- fit_mrs(mrs_data, basis, method = "abfit", time = FALSE)
  
  expected <- 7.313048e-5
  expect_equal(fit_res$res_tab$res.deviance, expected, tolerance = 1e-5,
               scale = expected)
})