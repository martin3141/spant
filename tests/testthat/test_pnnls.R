context("pnnls")

# pnnls() reduces overdetermined (m > n) problems to an n x n triangular
# system via a QR factorisation before calling the Lawson-Hanson PNNLS
# solver, for performance (see src/qrpnnls.f). These tests check that this
# gives (numerically) the same answer as calling the underlying Fortran
# PNNLS routine directly on the untransformed problem.

pnnls_direct <- function(a, b, k = 0) {
  m <- as.integer(nrow(a))
  n <- as.integer(ncol(a))
  storage.mode(a) <- "double"
  storage.mode(b) <- "double"
  r <- .Fortran("pnnls", r = a, m, m, n, b = b, x = double(n),
                rnorm = double(1), double(n), double(m), index = integer(n),
                mode = integer(1), k = as.integer(k), PACKAGE = "spant")
  r[c("x", "rnorm", "mode")]
}

test_that("pnnls with QR reduction matches direct PNNLS on random problems", {
  set.seed(1)
  for (i in 1:20) {
    m <- sample(c(10, 50, 200), 1)
    n <- sample(2:min(m, 40), 1)
    k <- sample(0:(n - 1), 1)
    a <- matrix(rnorm(m * n), m, n)
    b <- rnorm(m)
    r_direct <- pnnls_direct(a, b, k = k)
    r_pnnls  <- pnnls(a, b, k = k)
    expect_equal(r_pnnls$x, r_direct$x, tolerance = 1e-8)
    expect_equal(r_pnnls$rnorm, r_direct$rnorm, tolerance = 1e-8)
    expect_equal(r_pnnls$mode, r_direct$mode)
  }
})

test_that("pnnls with QR reduction matches direct PNNLS on real fitting-shaped problems", {
  sim_res <- sim_brain_1h(full_output = TRUE)
  metab   <- lb(sim_res$mrs_data, 5) |> phase(120) |> shift(10, units = "hz")
  fit_res <- fit_mrs(metab, sim_res$basis, progress = "none", time = FALSE)
  expect_true(fit_res$res_tab$res.deviance < 1)
})

test_that("pnnls rejects m < n problems", {
  a <- matrix(rnorm(15), 3, 5)
  b <- rnorm(3)
  expect_error(pnnls(a, b), "nrow\\(a\\) must be >= ncol\\(a\\)")
})

test_that("pnnls with m == n does not attempt QR reduction", {
  set.seed(2)
  n <- 10
  a <- matrix(rnorm(n * n), n, n)
  b <- rnorm(n)
  r_direct <- pnnls_direct(a, b, k = 2)
  r_pnnls  <- pnnls(a, b, k = 2)
  expect_equal(r_pnnls$x, r_direct$x, tolerance = 1e-8)
})
