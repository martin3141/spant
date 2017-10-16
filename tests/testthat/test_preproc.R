context("preprocessing")

test_that("first point scaling is correct", {
  expect_equal(sim_resonances()$data[1], 0.5 + 0i)
})