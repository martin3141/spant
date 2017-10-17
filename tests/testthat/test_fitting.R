context("fitting")

test_that("a full fitting pipeline works", {
  fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")
  mrs_data <- read_mrs(fname, format = "spar_sdat")
  mrs_proc <- hsvd_filt(mrs_data)
  mrs_proc <- align(mrs_proc, 2.01)
  basis <- sim_basis_1h_brain_press(mrs_proc)
  fit_res <- fit_mrs(mrs_proc, basis)
  # fit time won't be stable, so best to remove
  fit_res$proc_time <- NULL
  expect_equal_to_reference(fit_res, "fit_res.rds", tolerance = 1e-4)
})