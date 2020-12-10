context("nifti")

test_that("nifti MRS data can be written and read back from disk", {
  
  set.seed(1)
  
  # X = 3, Y = 5, Z = 7, dynamic = 16, coil = 8, FID = 256
  mrs_dims <- c(1, 3, 5, 7, 16, 8, 256)
  data <- rnorm(prod(mrs_dims)) + 1i * rnorm(prod(mrs_dims))
  dim(data) <- mrs_dims
  sim_mrs <- array2mrs_data(data)
  sim_mrs$resolution[2:7] <- c(9, 2, 3, 4, NA, 5e-4)
  
  # slice is rotated at an angle of 30 degrees in the z-plane
  sim_mrs$row_vec <- c(3 ^ 0.5 / 2, 1 / 2, 0)
  sim_mrs$col_vec <- c(1 / 2, -3 ^ 0.5 / 2, 0)
  sim_mrs$pos_vec <- c(0.5, -0.6, 1)
  
  sim_mrs$affine <- diag(1, 4, 4)
  
  tempf <- tempfile(fileext = ".nii.gz")
  write_mrs(sim_mrs, tempf)
  
  sim_mrs_nii <- read_mrs(tempf)
  
  # TODO add nucleus to checks
  expect_equal(sim_mrs$data, sim_mrs_nii$data)
  expect_equal(sim_mrs$resolution, sim_mrs_nii$resolution)
  expect_equal(sim_mrs$ft, sim_mrs_nii$ft)
  expect_equal(sim_mrs$ref, sim_mrs_nii$ref)
  expect_equal(sim_mrs$freq_domain, sim_mrs_nii$freq_domain)
  expect_equal(sim_mrs$te, sim_mrs_nii$te)
  #expect_equal(sim_mrs$row_vec, sim_mrs_nii$row_vec, tolerance = 1e-6)
  #expect_equal(sim_mrs$col_vec, sim_mrs_nii$col_vec, tolerance = 1e-6)
  #expect_equal(sim_mrs$pos_vec, sim_mrs_nii$pos_vec, tolerance = 1e-6)
  expect_equal(sim_mrs$affine, sim_mrs_nii$affine, tolerance = 1e-6)
})