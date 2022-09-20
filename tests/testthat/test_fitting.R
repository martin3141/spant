context("fitting")

test_that("Test data stimulation step", {
  
  # simulate some data and check consistency to a reference file
  sim_res <- sim_brain_1h(pul_seq = seq_press_ideal, full_output = TRUE)
  
  expect_equal_to_reference(sim_res, "abfit_sim_res.rds", tolerance = 1e-6)
  
  # create a realistic looking spectrum
  set.seed(1)
  
  mrs_data <- sim_res$mrs_data |> lb(6) |> scale_mrs_amp(1 / 10) |>
              phase(130) |> shift(0.1) |> add_noise(0.001, fd = FALSE)
  
  expect_equal_to_reference(mrs_data, "abfit_sim_mrs_data.rds",
                            tolerance = 1e-6)
})

test_that("Test ABfit with default options", {
  mrs_data <- readRDS("abfit_sim_mrs_data.rds")
  sim_res  <- readRDS("abfit_sim_res.rds")
  
  fit_res  <- fit_mrs(mrs_data, sim_res$basis, method = "abfit", time = FALSE,
                      progress = "none")
  
  expect_equal_to_reference(fit_res, "abfit_res_default.rds", tolerance = 1e-1)
  ref_val <- 0.4465198
  expect_equal(fit_res$res_tab$res.deviance, ref_val, tolerance = 1.5e-3)
})

test_that("Test ABfit without iterative optimisation", {
  mrs_data <- readRDS("abfit_sim_mrs_data.rds")
  sim_res  <- readRDS("abfit_sim_res.rds")
  
  # don't run any optimisation steps
  opts     <- abfit_opts(maxiters_pre = 0, maxiters = 0)
  
  # undo the large phase error to avoid FWHM measurement warnings
  fit_res  <- fit_mrs(mrs_data |> phase(-130), sim_res$basis, method = "abfit",
                      opts = opts, time = FALSE, progress = "none")
  
  expect_equal_to_reference(fit_res, "abfit_res_no_optim.rds", tolerance = 1e-5)
})

test_that("Test ABfit coarse-fitting steps only", {
  mrs_data <- readRDS("abfit_sim_mrs_data.rds")
  sim_res  <- readRDS("abfit_sim_res.rds")
  
  # don't run the fine fitting step
  opts <- abfit_opts(maxiters = 0)
  fit_res <- fit_mrs(mrs_data, sim_res$basis, method = "abfit", opts = opts,
                     time = FALSE, progress = "none")
  
  # this test is a bit more sensitive to differences between platforms
  expect_equal_to_reference(fit_res, "abfit_res_coarse.rds", tolerance = 2e-1)
  ref_val <- 0.2556892
  expect_equal(fit_res$res_tab$res.deviance, ref_val, tolerance = 1e-3)
})

test_that("Test ABfit fine-fitting only", {
  mrs_data <- readRDS("abfit_sim_mrs_data.rds")
  sim_res  <- readRDS("abfit_sim_res.rds")
  
  # don't run the pre fitting step
  opts <- abfit_opts(maxiters_pre = 0)
  
  # undo the large phase error to avoid FWHM measurement warnings
  fit_res  <- fit_mrs(mrs_data |> phase(-120), sim_res$basis, method = "abfit",
                      opts = opts, time = FALSE, progress = "none")
  
  expect_equal_to_reference(fit_res, "abfit_res_fine.rds", tolerance = 1e-6)
})

