beta2lw <- function(beta) {2 * (-beta * log(0.5)) ^ 0.5 / pi}

lw2beta <- function(lw) {(lw * pi / 2) ^ 2 / (-log(0.5))}

alpha2lw <- function(alpha) {alpha / pi}

lw2alpha <- function(lw) {lw * pi}

ft_shift <- function(vec_in) {pracma::fftshift(stats::fft(vec_in))}

ift_shift <- function(vec_in) {pracma::ifft(pracma::ifftshift(vec_in))}

hilbert <- function(x)
{
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