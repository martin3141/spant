basis <- sim_basis_1h_brain_press()

test_that("default basis is consistant", {
  expect_equal_to_reference(basis, "def_basis.rds", tolerance = 1e-6)
})

test_that("LCM .basis files can be written and read back", {
  basis_f <- tempfile()
  write_basis(basis, basis_f)
  basis_read <- read_basis(basis_f, sort_basis = FALSE)
  expect_equal_to_reference(basis_read, "def_basis.rds", tolerance = 1e-6)
})