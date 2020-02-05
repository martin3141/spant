context("fitting")

test_that("a full fitting pipeline works", {
  fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
                       package = "spant")
  
  mrs_data <- read_mrs(fname, format = "spar_sdat")
  mrs_proc <- hsvd_filt(mrs_data)
  mrs_proc <- align(mrs_proc, 2.01)
  basis <- sim_basis_1h_brain_press(mrs_proc)
  fit_res <- fit_mrs(mrs_proc, basis, method = "varpro_3p", time = FALSE)
  expect_equal_to_reference(fit_res, "fit_res_varpro_3p.rds", tolerance = 1e-4)
  
  fit_res <- fit_mrs(mrs_proc, basis, method = "varpro", time = FALSE)
  expect_equal_to_reference(fit_res, "fit_res_varpro.rds", tolerance = 1e-4)
})