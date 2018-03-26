context("qm simulator")

test_that("default basis is consistant", {
  expect_equal_to_reference(sim_basis_1h_brain_press(), "def_basis.rds",
                            tolerance = 1e-6)
})