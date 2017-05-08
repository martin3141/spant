beta2lw <- function(beta) {2 * (-beta * log(0.5)) ^ 0.5 / pi}

lw2beta <- function(lw) {(lw * pi / 2) ^ 2 / (-log(0.5))}

alpha2lw <- function(alpha) {alpha / pi}

lw2alpha <- function(lw) {lw * pi}

ft_shift <- function(vec_in) {pracma::fftshift(stats::fft(vec_in))}

ift_shift <- function(vec_in) {pracma::ifft(pracma::ifftshift(vec_in))}

hilbert <- function(x) {
  x <- Re(x)
  x <- stats::fft(x)
  n <- length(x)
  h <- vector()
  h[1] <- 1
  h[floor(n / 2) + 1] <- 1
  h[2:round(n / 2 + 0.5)] <- 2
  h[floor(n / 2 + 2):n] <- 0
  y <- h * x
  stats::fft(y, inverse = TRUE) / n
}

is.installed <- function(mypkg) {
  is.element(mypkg, utils::installed.packages()[,1]) 
}

rotate_vec <- function(vec_in, ax, theta) {
  ct <- cos(theta)
  st <- sin(theta)
  rotate_mat <- matrix(nrow = 3, ncol = 3)
  rotate_mat[1, 1] = ct + ax[1] * ax[1] * (1 - ct)
  rotate_mat[1, 2] = ax[1] * ax[2] * (1 - ct) - ax[3] * st
  rotate_mat[1, 3] = ax[1] * ax[3] * (1 - ct) + ax[2] * st
  rotate_mat[2, 1] = ax[2] * ax[1] * (1 - ct) + ax[3] * st
  rotate_mat[2, 2] = ct + ax[1] * ax[1] * (1 - ct)
  rotate_mat[2, 3] = ax[2] * ax[3] * (1 - ct) - ax[1] * st
  rotate_mat[3, 1] = ax[3] * ax[1] * (1 - ct) - ax[2] * st
  rotate_mat[3, 2] = ax[3] * ax[2] * (1 - ct) + ax[1] * st
  rotate_mat[3, 3] = ct + ax[3] * ax[3] * (1 - ct)
  vec_out <- rotate_mat %*% vec_in
  as.numeric(vec_out / sum((vec_out) ^ 2))
}

# From : http://mathworld.wolfram.com/LeastSquaresFittingExponential.html
#' @export
measure_lorentz_amp <- function(y, t, start_pnt = 10, end_pnt = 50) {
  # crop to time region
  y <- y[start_pnt:end_pnt]
  t <- t[start_pnt:end_pnt]
  
  a <- sum((t ^ 2) * y) * sum(y * log(y)) -  sum(t * y) * sum(t * y * log(y))
  a <- a / (sum(y) * sum(t ^ 2 * y) - sum(t * y) ^ 2)
  A <- exp(a)
  A
}