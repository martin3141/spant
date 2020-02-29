context("fitting")

test_that("a full fitting pipeline works", {
  fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
                       package = "spant")
  
  mrs_data <- read_mrs(fname, format = "spar_sdat")
  basis    <- sim_basis_1h_brain_press(mrs_data)
  
  fit_res  <- fit_mrs(mrs_data, basis, method = "abfit", time = FALSE)
  expect_equal_to_reference(fit_res, "fit_res_abfit.rds", tolerance = 1e-4)
})